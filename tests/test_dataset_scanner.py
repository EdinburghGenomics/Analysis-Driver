import os
from unittest.mock import Mock, patch
from egcg_core.constants import DATASET_NEW, DATASET_ABORTED
from analysis_driver.dataset import RunDataset, SampleDataset
from analysis_driver.dataset_scanner import DatasetScanner, RunScanner, SampleScanner
from tests.test_analysisdriver import TestAnalysisDriver
from tests.test_dataset import patched_expected_yield, fake_proc, seed_directories, clean


class FakeEntity(Mock):
    def __init__(self, name, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.name = name


def ppath(*parts):
    return '.'.join(('analysis_driver', 'dataset_scanner') + parts)


def patched_get(content=None):
    if content is None:
        content = [fake_proc]
    return patch(ppath('rest_communication', 'get_documents'), return_value=content)


class TestScanner(TestAnalysisDriver):
    def setUp(self):
        self.base_dir = os.path.join(self.assets_path, 'dataset_scanner')
        seed_directories(self.base_dir)
        self._setup_scanner()
        assert self.scanner.input_dir == self.base_dir

        self.triggerignore = os.path.join(self.scanner.input_dir, '.triggerignore')
        with open(self.triggerignore, 'w') as f:
            f.write('ignored_dataset\n')

    def tearDown(self):
        clean(self.base_dir)

    def test_get_datasets_for_status(self, *args):
        pass

    def test_report(self):
        dsets = {'new': ['this'], 'ready': ['that', 'other'], 'failed': ['another']}
        captured_stdout = []
        patched_scan = patch(ppath('DatasetScanner.scan_datasets'), return_value=dsets)
        with patched_scan, patch('builtins.print', new=captured_stdout.append):
            self.scanner.report(all_datasets=True)

        expected = [
            '========= ' + self.scanner.__class__.__name__ + ' report =========',
            'dataset location: ' + self.base_dir,
            '=== new ===',
            'this',
            '=== ready ===',
            'that',
            'other',
            '=== failed ===',
            'another',
            '__________________________________________'
        ]
        observed = captured_stdout[0].split('\n')
        print(observed)
        print(expected)
        assert observed == expected

    def test_triggerignore(self):
        with open(self.triggerignore, 'r') as f:
            assert f.readlines() == ['ignored_dataset\n']
            assert os.path.isdir(os.path.join(self.base_dir, 'ignored_dataset'))
            assert self.scanner._triggerignore == ['ignored_dataset']

    def _setup_scanner(self):
        self.scanner = DatasetScanner({'input_dir': self.base_dir})

    def _assert_datasets_equal(self, obs, exp):
        assert obs.name == exp.name
        assert obs.path == exp.path
        assert obs.use_int_dir == exp.use_int_dir


class TestRunScanner(TestScanner):
    def test_get_dataset(self):
        test_dataset_path = os.path.join(self.base_dir, 'test_dataset')
        os.mkdir(test_dataset_path)
        with patched_get():
            obs = self.scanner.get_dataset('test_dataset')
            exp = RunDataset('test_dataset', test_dataset_path, self.scanner.use_int_dir)
        self._assert_datasets_equal(obs, exp)

    def test_get_dataset_records_for_status(self):
        for d in self.scanner.expected_bcl_subdirs:
            os.makedirs(os.path.join(self.base_dir, 'test_dataset', d))
        with patched_get([{'run_id': 'test_dataset'}]):
            obs = self.scanner._get_dataset_records_for_status(DATASET_NEW)
        assert obs[0]['run_id'] == 'test_dataset'

        fake_data = {
            'dataset_type': 'run',
            'dataset_name': 'test_dataset',
            'status': DATASET_ABORTED,
            'run_id': 'test_dataset'
        }
        with patched_get([fake_data]):
            obs = self.scanner._get_dataset_records_for_status(DATASET_ABORTED)[0]['run_id']
            assert obs == 'test_dataset'

    def test_datasets_on_disk(self):
        for d in self.scanner.expected_bcl_subdirs:
            os.makedirs(os.path.join(self.base_dir, 'test_dataset', d))
        obs = self.scanner._datasets_on_disk()
        exp = ['test_dataset']
        assert obs == exp

    def test_scan_datasets(self):
        fake_datasets = []
        with patched_get():
            for x in ('this', 'that', 'other'):
                d = RunDataset(x, os.path.join(self.base_dir, x), False)
                assert d.most_recent_proc.entity
                fake_datasets.append(d)
        with patch(ppath('DatasetScanner._get_datasets_for_status'), return_value=fake_datasets) as p:
            obs = self.scanner.scan_datasets(DATASET_NEW)
            p.assert_any_call(DATASET_NEW)
            exp = {DATASET_NEW: '[other, that, this]'}
            assert str(obs[DATASET_NEW]) == exp[DATASET_NEW]

    def _setup_scanner(self):
        self.scanner = RunScanner({'input_dir': self.base_dir})


class TestSampleScanner(TestScanner):
    def _setup_scanner(self):
        self.scanner = SampleScanner({'input_dir': self.base_dir})

    def test_get_dataset(self):
        with patched_expected_yield(), patched_get() as p:
            observed = self.scanner.get_dataset('test_dataset')
            expected = SampleDataset('test_dataset')
            assert observed.name == expected.name
            assert observed.data_threshold == expected.data_threshold
            assert observed.run_elements == [fake_proc]
            p.assert_called_with('run_elements', where={'sample_id': 'test_dataset', 'useable': 'yes'})

    @patched_get([{'sample_id': 'a_sample_id'}])
    def test_get_dataset_records_for_status(self, mocked_get):
        assert self.scanner._get_dataset_records_for_status(DATASET_NEW)[0] == {'sample_id': 'a_sample_id'}
        mocked_get.assert_any_call(
            'aggregate/samples',
            match={'proc_status': None},
            paginate=False
        )

    @patched_get([{'sample_id': 'a_sample_id'}])
    def test_get_datasets_for_status(self, mocked_get):
        patched_list_samples = patch(
            ppath('get_list_of_samples'),
            return_value=[FakeEntity('a_sample_id', udf={'Yield for Quoted Coverage (Gb)': 1})]
        )
        with patched_list_samples:
            d = self.scanner._get_datasets_for_status(DATASET_NEW)[0]
        assert d.name == 'a_sample_id'
        assert d._data_threshold == 1000000000
        mocked_get.assert_any_call('aggregate/samples', match={'proc_status': None}, paginate=False)
        assert d.run_elements == [{'sample_id': 'a_sample_id'}]
        mocked_get.assert_any_call('run_elements', where={'useable': 'yes', 'sample_id': 'a_sample_id'})

    def test_scan_datasets(self):
        fake_datasets = {DATASET_NEW: []}
        with patched_get():
            for x in ('other', 'that', 'this'):
                d = SampleDataset(x)
                assert d.most_recent_proc.entity
                fake_datasets[DATASET_NEW].append(d)

        def fake_get_datasets(*args):
            return fake_datasets[args[1]]

        patched_is_ready = patch(ppath('SampleDataset._is_ready'), return_value=False)
        patched_get_datasets = patch(ppath('SampleScanner._get_datasets_for_status'), new=fake_get_datasets)
        with patched_is_ready, patched_get_datasets:
            assert self.scanner.scan_datasets(DATASET_NEW) == fake_datasets
