import os
import shutil
import pytest
from unittest.mock import patch, PropertyMock
from tests.test_analysisdriver import TestAnalysisDriver
from tests.test_clarity import FakeEntity
from analysis_driver.exceptions import AnalysisDriverError
from analysis_driver.dataset_scanner import DatasetScanner, RunScanner, SampleScanner, Dataset, RunDataset,\
    SampleDataset, MostRecentProc
from analysis_driver.constants import DATASET_NEW, DATASET_READY, DATASET_PROCESSING, DATASET_PROCESSED_FAIL,\
    DATASET_PROCESSED_SUCCESS, DATASET_ABORTED, DATASET_REPROCESS, DATASET_FORCE_READY

directories_to_create = ('dataset_ready', 'dataset_not_ready', 'ignored_dataset')
ready_datasets = ('dataset_ready',)


def seed_directories(base_dir):
    for d in directories_to_create:
        os.makedirs(os.path.join(base_dir, d), exist_ok=True)
        if d in ready_datasets:
            _touch(os.path.join(base_dir, d, 'RTAComplete.txt'))


def clean(base_dir):
    shutil.rmtree(os.path.join(base_dir))


def _touch(path):
    open(path, 'w').close()


def ppath(*parts):
    return '.'.join(('analysis_driver', 'dataset_scanner') + parts)


patched_patch = patch(ppath('rest_communication', 'patch_entry'))
patched_post = patch(ppath('rest_communication', 'post_entry'))
patched_pid = patch(ppath('os.getpid'), return_value=1)
patched_update = patch(ppath('MostRecentProc.update_entity'))
patched_initialise = patch(ppath('MostRecentProc.initialise_entity'))
patched_finish = patch(ppath('MostRecentProc.finish'))

patched_ntf_start = patch(ppath('ntf.start_pipeline'))
patched_ntf_end = patch(ppath('ntf.end_pipeline'))
patched_ntf_start_stage = patch(ppath('ntf.start_stage'))
patched_ntf_end_stage = patch(ppath('ntf.end_stage'))


def patched_expected_yield(y=1000000000):
    return patch(ppath('get_expected_yield_for_sample'), return_value=y)

patched_stages = patch(
    ppath('Dataset.stages'),
    new_callable=PropertyMock(return_value=['this', 'that', 'other'])
)

fake_proc = {
    'proc_id': 'test_test_now',
    'dataset_type': 'test',
    'dataset_name': 'test'
}


def patched_get(content=None):
    if content is None:
        content = [fake_proc]
    return patch(ppath('rest_communication', 'get_documents'), return_value=content)


def patched_datetime(time='now'):
    return patch(ppath('MostRecentProc._now'), return_value=time)


class TestMostRecentProc(TestAnalysisDriver):
    def setUp(self):
        with patched_get():
            self.proc = MostRecentProc('test', 'test')

    def init_proc_with_fake_data(self):
        self.proc._rest_entity = fake_proc.copy()
        self.proc.local_entity = fake_proc.copy()

    def test_rest_entity_pre_existing(self):
        self.init_proc_with_fake_data()
        assert self.proc.rest_entity == fake_proc

    @patched_get()
    @patched_initialise
    def test_rest_entity_not_pre_existing(self, mocked_initialise, mocked_get):
        with patched_datetime():
            self.proc._rest_entity = None
            x = self.proc.rest_entity
            assert x == fake_proc
            mocked_get.assert_called_with(
                'analysis_driver_procs',
                where={'dataset_type': 'test', 'dataset_name': 'test'},
                sort='-_created'
            )

        self.proc._rest_entity = None
        with patched_get([]):
            y = self.proc.rest_entity
            assert y is None
            mocked_get.assert_called_with(
                'analysis_driver_procs',
                where={'dataset_type': 'test', 'dataset_name': 'test'},
                sort='-_created'
            )
            assert mocked_initialise.call_count == 1

    @patched_post
    @patched_patch
    def test_initialise_entity(self, mocked_patch, mocked_post):
        assert self.proc._rest_entity == fake_proc
        with patched_datetime():
            self.proc.initialise_entity()
        mocked_post.assert_called_with(
            'analysis_driver_procs',
            {'proc_id': 'test_test_now', 'dataset_type': 'test', 'dataset_name': 'test'}
        )
        mocked_patch.assertcalled_with(
            'tests',
            {'analysis_driver_procs': ['test_test_now']},
            id_field='test_id',
            element_id='test',
            update_lists=['analysis_driver_procs']
        )

        with patched_datetime('then'):
            self.proc.initialise_entity()
        mocked_post.assert_called_with(
            'analysis_driver_procs',
            {'proc_id': 'test_test_then', 'dataset_type': 'test', 'dataset_name': 'test'}
        )
        mocked_patch.assert_called_with(
            'tests',
            {'analysis_driver_procs': ['test_test_then']},
            id_field='test_id',
            element_id='test',
            update_lists=['analysis_driver_procs']
        )

    @patched_patch
    def test_sync(self, mocked_patch):
        self.init_proc_with_fake_data()
        self.proc.sync()
        assert self.proc._rest_entity == {
            'proc_id': 'test_test_now',
            'dataset_type': 'test',
            'dataset_name': 'test'
        }
        assert self.proc.local_entity == self.proc.rest_entity

        self.proc.local_entity.update({'this': 'that'})
        self.proc.sync()
        mocked_patch.assert_called_with(
            'analysis_driver_procs',
            {'this': 'that'},
            id_field='proc_id',
            element_id='test_test_now'
        )
        assert self.proc._rest_entity == {
            'proc_id': 'test_test_now',
            'dataset_type': 'test',
            'dataset_name': 'test',
            'this': 'that'
        }

    @patch(ppath('MostRecentProc.sync'))
    def test_update_entity(self, mocked_sync):
        self.init_proc_with_fake_data()
        self.proc.update_entity(other='another')
        assert self.proc.local_entity == {
            'proc_id': 'test_test_now',
            'dataset_type': 'test',
            'dataset_name': 'test',
            'other': 'another'
        }
        mocked_sync.assert_called()

    @patched_update
    def test_start(self, mocked_update):
        self.init_proc_with_fake_data()
        with patched_pid:
            self.proc.start()
        mocked_update.assert_called_with(status=DATASET_PROCESSING, pid=1)

    @patched_update
    def test_finish(self, mocked_update):
        self.init_proc_with_fake_data()
        with patched_datetime():
            self.proc.finish(DATASET_PROCESSED_SUCCESS)
        mocked_update.assert_called_with(status=DATASET_PROCESSED_SUCCESS, pid=None, end_date='now')

    @patched_update
    def test_start_stage(self, mocked_update):
        self.init_proc_with_fake_data()
        with patched_datetime('then'):
            self.proc.start_stage('test_stage')
        mocked_update.assert_called_with(
            stages=[
                {'date_started': 'then', 'stage_name': 'test_stage'}
            ]
        )

        self.proc.local_entity['stages'] = [{'date_started': 'then', 'stage_name': 'test_stage'}]
        with patched_datetime('later'):
            self.proc.start_stage('another_stage')
        mocked_update.assert_called_with(
            stages=[
                {'date_started': 'then', 'stage_name': 'test_stage'},
                {'date_started': 'later', 'stage_name': 'another_stage'}
            ]
        )

    @patched_update
    def test_end_stage(self, mocked_update):
        self.init_proc_with_fake_data()
        with pytest.raises(KeyError) as e:
            self.proc.end_stage('test_stage')
            assert str(e) == 'stages'

        self.proc.local_entity['stages'] = [
            {'date_started': 'now', 'stage_name': 'this'},
            {'date_started': 'then', 'stage_name': 'this'},
            {'date_started': 'later', 'stage_name': 'that'}
        ]
        with patched_datetime('finally'):
            self.proc.end_stage('this')
        mocked_update.assert_called_with(
            stages=[
                {'date_started': 'now', 'stage_name': 'this', 'date_finished': 'finally', 'exit_status': 0},
                {'date_started': 'then', 'stage_name': 'this', 'date_finished': 'finally', 'exit_status': 0},
                {'date_started': 'later', 'stage_name': 'that'}
            ]
        )


class TestDataset(TestAnalysisDriver):
    def setUp(self):
        self.base_dir = os.path.join(self.assets_path, 'dataset_scanner')
        seed_directories(self.base_dir)
        self.setup_dataset()

    def tearDown(self):
        clean(self.base_dir)

    def test_dataset_status(self):
        assert self.dataset.most_recent_proc.local_entity.get('status') is None
        assert self.dataset.dataset_status == DATASET_NEW

        self.dataset.most_recent_proc.local_entity['status'] = 'a_status'
        assert self.dataset.dataset_status == 'a_status'

    def test_stages(self):
        assert self.dataset.stages == []
        self.dataset.most_recent_proc.local_entity['stages'] = [
            {'stage_name': 'a_stage', 'date_started': 'now', 'date_finished': 'finally'},
            {'stage_name': 'another_stage', 'date_started': 'then'},
            {'stage_name': 'yet_another_stage', 'date_started': 'finally'}
        ]
        assert self.dataset.stages == ['another_stage', 'yet_another_stage']

    @patched_ntf_start
    @patched_initialise
    @patch(ppath('MostRecentProc.start'))
    def test_start(self, mocked_start, mocked_update, mocked_ntf):
        self.dataset.most_recent_proc.local_entity['status'] = DATASET_PROCESSING
        with pytest.raises(AssertionError):
            self.dataset.start()

        del self.dataset.most_recent_proc.local_entity['status']
        self.dataset.start()
        for m in (mocked_update, mocked_ntf, mocked_start):
            assert m.call_count == 1

    @patched_ntf_end
    @patched_finish
    def test_succeed(self, mocked_finish, mocked_ntf):
        with pytest.raises(AssertionError):
            self.dataset.succeed()

        self.dataset.most_recent_proc.local_entity['status'] = DATASET_PROCESSING
        self.dataset.succeed()
        mocked_ntf.assert_called_with(0)
        mocked_finish.assert_called_with(DATASET_PROCESSED_SUCCESS)

    @patched_ntf_end
    @patched_finish
    def test_fail(self, mocked_finish, mocked_ntf):
        with pytest.raises(AssertionError):
            self.dataset.succeed()

        self.dataset.most_recent_proc.local_entity['status'] = DATASET_PROCESSING
        self.dataset.fail(1)
        mocked_ntf.assert_called_with(1)
        mocked_finish.assert_called_with(DATASET_PROCESSED_FAIL)

    @patched_finish
    def test_abort(self, mocked_finish):
        self.dataset.abort()
        mocked_finish.assert_called_with(DATASET_ABORTED)

    @patch(ppath('MostRecentProc.change_status'))
    def test_reset(self, mocked_change_status):
        self.dataset.reset()
        mocked_change_status.assert_called_with(DATASET_REPROCESS)

    @patched_ntf_start_stage
    @patch(ppath('MostRecentProc.start_stage'))
    def test_start_stage(self, mocked_start_stage, mocked_ntf):
        self.dataset.start_stage('a_stage')
        mocked_start_stage.assert_called_with('a_stage')
        mocked_ntf.assert_called_with('a_stage')

    @patched_ntf_end_stage
    @patch(ppath('MostRecentProc.end_stage'))
    def test_end_stage(self, mocked_end_stage, mocked_ntf):
        self.dataset.end_stage('a_stage')
        mocked_end_stage.assert_called_with('a_stage', 0)
        mocked_ntf.assert_called_with('a_stage', 0)

    def test_str(self):
        with patched_stages:
            assert str(self.dataset) == 'test_dataset -- this, that, other'

    def setup_dataset(self):
        with patched_get():
            self.dataset = _TestDataset(
                'test_dataset',
                MostRecentProc(
                    'test',
                    'test',
                    initial_content={'date_started': 'now', 'dataset_name': 'None', 'dataset_type': 'None'}
                )
            )


class _TestDataset(Dataset):
    type = 'None'
    endpoint = 'None'
    id_field = 'None'

    def _is_ready(self):
        pass


class TestRunDataset(TestDataset):
    def test_is_ready(self):
        datasets = (
            ('dataset_ready', True),
            ('dataset_not_ready', False)
        )
        with patched_get():
            for d_name, rta_complete in datasets:
                d = RunDataset(d_name, os.path.join(self.base_dir, d_name), use_int_dir=False)
                assert d._is_ready() == rta_complete

    def test_dataset_status(self):
        super().test_dataset_status()
        del self.dataset.most_recent_proc.local_entity['status']
        assert not self.dataset.rta_complete()
        assert self.dataset.dataset_status == DATASET_NEW
        os.mkdir(os.path.join(self.base_dir, self.dataset.name))
        _touch(os.path.join(self.base_dir, self.dataset.name, 'RTAComplete.txt'))
        assert self.dataset.rta_complete()
        assert self.dataset.dataset_status == DATASET_READY

    def setup_dataset(self):
        self.dataset = RunDataset(
            'test_dataset',
            os.path.join(self.base_dir, 'test_dataset'),
            use_int_dir=False,
            most_recent_proc=MostRecentProc(
                'test',
                'test',
                initial_content={'date_started': 'now', 'dataset_name': 'None', 'dataset_type': 'None'}
            )
        )


class TestSampleDataset(TestDataset):
    def test_dataset_status(self):
        with patched_expected_yield():
            super().test_dataset_status()

    @patch(ppath('MostRecentProc.change_status'))
    def test_force(self, mocked_change_status):
        self.dataset.force()
        mocked_change_status.assert_called_with(DATASET_FORCE_READY)

    @patched_get()
    def test_read_data(self, mocked_get):
        assert self.dataset._read_data() == [fake_proc]
        mocked_get.assert_called_with('run_elements', where={'sample_id': 'test_dataset', 'useable': 'yes'})

    def test_amount_data(self):
        assert self.dataset._amount_data() == 480

    def test_runs(self):
        with patched_get(self.dataset.run_elements):
            assert self.dataset._runs() == ['a_run_id', 'another_run_id']

    @patched_expected_yield()
    def test_data_threshold(self, mocked_exp_yield):
        self.dataset._data_threshold = None
        assert self.dataset.data_threshold == 1000000000
        mocked_exp_yield.assert_called_with('test_dataset')

    @patched_expected_yield(None)
    def test_no_data_threshold(self, mocked_exp_yield):
        self.dataset._data_threshold = None
        with pytest.raises(AnalysisDriverError) as e:
            _ = self.dataset.data_threshold
        assert 'Could not find data threshold in LIMS' in str(e)
        mocked_exp_yield.assert_called_with('test_dataset')

    @patched_expected_yield()
    def test_is_ready(self, mocked_instance):
        self.dataset._data_threshold = None
        assert not self.dataset._is_ready()
        self.dataset._run_elements = [
            {
                'clean_q30_bases_r1': 1200000000,
                'clean_q30_bases_r2': 1150000000
            },
            {
                'clean_q30_bases_r1': 1100000000,
                'clean_q30_bases_r2': 1350000000
            }
        ]
        assert self.dataset._is_ready()
        assert mocked_instance.call_count == 1  # even after 2 calls to data_threshold

    def test_str(self):
        expected_str = 'test_dataset -- this, that, other  (480 / 1000000000  from a_run_id, another_run_id) '
        self.dataset._data_threshold = None
        with patched_get(self.dataset.run_elements), patched_expected_yield(), patched_stages:
            print(expected_str)
            print(str(self.dataset))
            assert str(self.dataset) == expected_str

    def setup_dataset(self):
        with patched_get():
            self.dataset = SampleDataset(
                'test_dataset',
                most_recent_proc=MostRecentProc(
                    'test',
                    'test',
                    initial_content={'date_started': 'now', 'dataset_name': 'None', 'dataset_type': 'None'}
                )
            )
        self.dataset._run_elements = [
            {
                'run_id': 'a_run_id',
                'clean_q30_bases_r1': 120,
                'clean_q30_bases_r2': 115
            },
            {
                'run_id': 'another_run_id',
                'clean_q30_bases_r1': 110,
                'clean_q30_bases_r2': 135
            }
        ]
        self.dataset._data_threshold = 1000000000


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
            '=== failed ===',
            'another',
            '=== new ===',
            'this',
            '=== ready ===',
            'that',
            'other',
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
                fake_datasets.append(RunDataset(x, os.path.join(self.base_dir, x), False))
        with patch(ppath('DatasetScanner._get_datasets_for_status'), return_value=fake_datasets) as p:
            obs = self.scanner.scan_datasets(DATASET_NEW)
            p.assert_any_call(DATASET_NEW)
            exp = {DATASET_NEW: '[other, that, this]'}
            assert str(obs[DATASET_NEW]) == exp[DATASET_NEW]

            obs = self.scanner.scan_datasets(DATASET_NEW, flatten=True)
            assert str(obs) == '[other, that, this]'

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
            match={'proc_status': DATASET_NEW}
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
        mocked_get.assert_any_call('aggregate/samples', match={'proc_status': 'new'})
        assert d.run_elements == [{'sample_id': 'a_sample_id'}]
        mocked_get.assert_any_call('run_elements', where={'useable': 'yes', 'sample_id': 'a_sample_id'})

    def test_scan_datasets(self):
        fake_datasets = {DATASET_NEW: []}
        with patched_get():
            for x in ('other', 'that', 'this'):
                fake_datasets[DATASET_NEW].append(SampleDataset(x))

        def fake_get_datasets(*args):
            return fake_datasets[args[1]]

        patched_is_ready = patch(ppath('SampleDataset._is_ready'), return_value=False)
        patched_get_datasets = patch(ppath('SampleScanner._get_datasets_for_status'), new=fake_get_datasets)
        with patched_is_ready, patched_get_datasets:
            assert self.scanner.scan_datasets(DATASET_NEW) == fake_datasets
            assert self.scanner.scan_datasets(DATASET_NEW, flatten=True) == fake_datasets[DATASET_NEW]
