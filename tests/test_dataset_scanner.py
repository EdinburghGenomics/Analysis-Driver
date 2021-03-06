import os
from unittest.mock import Mock, patch
from egcg_core.constants import DATASET_READY, DATASET_NEW, DATASET_ABORTED
from analysis_driver.dataset import RunDataset, SampleDataset
from analysis_driver.dataset_scanner import DatasetScanner, RunScanner, SampleScanner
from tests.test_analysisdriver import TestAnalysisDriver, NamedMock
from tests.test_dataset import fake_proc


class FakeRunDataset(RunDataset):
    @property
    def running_stages(self):
        return []


def ppath(*parts):
    return '.'.join(('analysis_driver', 'dataset_scanner') + parts)


def patched_get(content=None):
    if content is None:
        content = [fake_proc]
    return patch(ppath('get_documents'), return_value=content)


class TestScanner(TestAnalysisDriver):
    def setUp(self):
        self.base_dir = os.path.join(self.assets_path, 'dataset_scanner')
        os.makedirs(self.base_dir, exist_ok=True)
        self._setup_scanner()
        self.scanner.input_dir = self.base_dir  # patch the config

    def test_report(self):
        def fake_dataset(name):
            return NamedMock(real_name=name, entity={'pid': 1}, report=Mock(return_value=name))

        dsets = {
            'new': [fake_dataset('this')],
            'ready': [fake_dataset('that'), fake_dataset('other')],
            'failed': [fake_dataset('another')]
        }
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
        assert observed == expected

    def test_triggerignore(self):
        triggerignore = os.path.join(self.scanner.input_dir, '.triggerignore')
        with open(triggerignore, 'w') as f:
            f.write('# a comment\nignored_dataset\n')

        with open(triggerignore, 'r') as f:
            assert f.readlines() == ['# a comment\n', 'ignored_dataset\n']
            assert self.scanner._triggerignore == ['ignored_dataset']

    def _setup_scanner(self):
        self.scanner = DatasetScanner()


class TestRunScanner(TestScanner):
    def test_get_dataset(self):
        obs = self.scanner.get_dataset('test_dataset')
        exp = FakeRunDataset('test_dataset')
        assert obs.name == exp.name

    def test_get_dataset_records_for_statuses(self):
        for d in self.scanner.expected_bcl_subdirs:
            os.makedirs(os.path.join(self.base_dir, 'test_dataset', d), exist_ok=True)
        with patched_get([{'run_id': 'test_dataset'}]):
            obs = self.scanner._get_dataset_records_for_statuses([DATASET_NEW])
        assert obs[0]['run_id'] == 'test_dataset'

        fake_data = {
            'dataset_type': 'run',
            'dataset_name': 'test_dataset',
            'status': DATASET_ABORTED,
            'run_id': 'test_dataset'
        }
        with patched_get([fake_data]):
            obs = self.scanner._get_dataset_records_for_statuses([DATASET_ABORTED])[0]['run_id']
            assert obs == 'test_dataset'

    def test_datasets_on_disk(self):
        for d in self.scanner.expected_bcl_subdirs:
            os.makedirs(os.path.join(self.base_dir, 'test_dataset', d), exist_ok=True)
        obs = self.scanner._datasets_on_disk()
        exp = ['test_dataset']
        assert obs == exp

    def test_scan_datasets(self):
        fake_datasets = []
        with patched_get(), patch('analysis_driver.dataset.rest_communication.get_documents'):
            for x in ('this', 'that', 'other'):
                d = FakeRunDataset(x, {'this': 'that'})
                assert d.most_recent_proc.entity
                fake_datasets.append(d)
        with patch(ppath('DatasetScanner._get_datasets_for_statuses'), return_value=fake_datasets) as p:
            obs = self.scanner.scan_datasets(DATASET_NEW)
            p.assert_any_call((DATASET_NEW, ))
            exp = {DATASET_READY: [FakeRunDataset(n) for n in ['other', 'that', 'this']]}
            assert [d.name for d in obs[DATASET_READY]] == [d.name for d in exp[DATASET_READY]]

    def _setup_scanner(self):
        self.scanner = RunScanner()


class TestSampleScanner(TestScanner):
    def _setup_scanner(self):
        self.scanner = SampleScanner()

    def test_get_dataset(self):
        observed = self.scanner.get_dataset('test_dataset')
        expected = SampleDataset('test_dataset')
        assert observed.name == expected.name

    @patched_get([{'sample_id': 'a_sample_id'}])
    def test_get_dataset_records_for_statuses(self, mocked_get):
        assert self.scanner._get_dataset_records_for_statuses([DATASET_NEW]) == [{'sample_id': 'a_sample_id'}]
        mocked_get.assert_any_call(
            'samples',
            where={'aggregated.most_recent_proc.status': None},
            all_pages=True,
            quiet=True,
            max_results=100
        )

    @patched_get([{'sample_id': 'a_sample_id'}])
    def test_get_datasets_for_statuses(self, mocked_get):
        d = self.scanner._get_datasets_for_statuses([DATASET_NEW])[0]
        assert d.name == 'a_sample_id'
        mocked_get.assert_any_call(
            'samples',
            where={'aggregated.most_recent_proc.status': None},
            all_pages=True,
            quiet=True,
            max_results=100
        )

    def test_scan_datasets(self):
        fake_datasets = {DATASET_NEW: []}
        for x in ('other', 'that', 'this'):
            s = SampleDataset(x, fake_proc)
            assert s.most_recent_proc.entity
            fake_datasets[DATASET_NEW].append(s)

        def fake_get_datasets(*args):
            return [d for arg in args[1] for d in fake_datasets.get(arg)]

        patched_is_ready = patch(ppath('SampleDataset._is_ready'), return_value=False)
        patched_get_datasets = patch(ppath('SampleScanner._get_datasets_for_statuses'), new=fake_get_datasets)
        with patched_is_ready, patched_get_datasets:
            assert self.scanner.scan_datasets(DATASET_NEW) == fake_datasets
