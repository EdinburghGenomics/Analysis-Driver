import os
import pytest
from shutil import rmtree
from unittest.mock import patch, Mock, PropertyMock
from tests.test_analysisdriver import TestAnalysisDriver
from egcg_core import constants as c
from analysis_driver.exceptions import AnalysisDriverError
from analysis_driver.dataset import Dataset, RunDataset, SampleDataset, MostRecentProc


def seed_directories(base_dir):
    for d in ('dataset_ready', 'dataset_not_ready', 'ignored_dataset'):
        os.makedirs(os.path.join(base_dir, d), exist_ok=True)
        if d in ('dataset_ready',):
            touch(os.path.join(base_dir, d, 'RTAComplete.txt'))


def clean(base_dir):
    rmtree(os.path.join(base_dir))


def touch(path):
    open(path, 'w').close()


def ppath(*parts):
    return '.'.join(('analysis_driver', 'dataset') + parts)


fake_proc = {'proc_id': 'test_test_now', 'dataset_type': 'test', 'dataset_name': 'test'}

patched_patch = patch(ppath('rest_communication', 'patch_entry'))
patched_post = patch(ppath('rest_communication', 'post_entry'))
patched_pid = patch(ppath('os.getpid'), return_value=1)
patched_update = patch(ppath('MostRecentProc.update_entity'))
patched_initialise = patch(ppath('MostRecentProc.initialise_entity'))
patched_finish = patch(ppath('MostRecentProc.finish'))
patched_stages = patch(
    ppath('Dataset.stages'),
    new_callable=PropertyMock(return_value=['this', 'that', 'other'])
)


def patched_get(content=None):
    if content is None:
        content = [fake_proc]
    return patch(ppath('rest_communication', 'get_documents'), return_value=content)


def patched_datetime(time='now'):
    return patch(ppath('MostRecentProc._now'), return_value=time)


def patched_expected_yield(y=1000000000):
    return patch(ppath('get_expected_yield_for_sample'), return_value=y)


def patched_stage_id(stage_id):
    return patch(ppath('MostRecentProc._stage_id'), return_value=stage_id)


class TestDataset(TestAnalysisDriver):
    def setUp(self):
        self.base_dir = os.path.join(self.assets_path, 'dataset_scanner')
        seed_directories(self.base_dir)
        self.setup_dataset()
        self.dataset.ntf = Mock()

    def tearDown(self):
        clean(self.base_dir)

    def test_dataset_status(self):
        assert self.dataset.most_recent_proc.entity.get('status') is None
        assert self.dataset.dataset_status == c.DATASET_NEW

        self.dataset.most_recent_proc.entity['status'] = 'a_status'
        assert self.dataset.dataset_status == 'a_status'

    def test_stages(self):
        assert self.dataset.stages == []
        self.dataset.most_recent_proc.entity['stages'] = [
            {'stage_name': 'a_stage', 'date_started': 'now', 'date_finished': 'finally'},
            {'stage_name': 'another_stage', 'date_started': 'then'},
            {'stage_name': 'yet_another_stage', 'date_started': 'finally'}
        ]
        assert self.dataset.stages == ['another_stage', 'yet_another_stage']

    @patched_initialise
    @patch(ppath('MostRecentProc.start'))
    @patch(ppath('MostRecentProc.retrieve_entity'))
    def test_start(self, mocked_retrieve_entity, mocked_start, mocked_update):
        self.dataset.most_recent_proc.entity['status'] = c.DATASET_PROCESSING
        with pytest.raises(AssertionError):
            self.dataset.start()

        del self.dataset.most_recent_proc.entity['status']
        self.dataset.start()
        for m in (mocked_update, self.dataset.ntf.start_pipeline, mocked_start):
            assert m.call_count == 1

    @patched_finish
    @patch(ppath('MostRecentProc.retrieve_entity'))
    def test_succeed(self, mocked_retrieve_entity, mocked_finish):
        with pytest.raises(AssertionError):
            self.dataset.succeed()

        self.dataset.most_recent_proc.entity['status'] = c.DATASET_PROCESSING
        self.dataset.succeed()
        self.dataset.ntf.end_pipeline.assert_called_with(0)
        mocked_finish.assert_called_with(c.DATASET_PROCESSED_SUCCESS)

    @patched_finish
    @patch(ppath('MostRecentProc.retrieve_entity'))
    def test_fail(self, mocked_retrieve_entity, mocked_finish):
        with pytest.raises(AssertionError):
            self.dataset.succeed()

        self.dataset.most_recent_proc.entity['status'] = c.DATASET_PROCESSING
        self.dataset.fail(1)
        self.dataset.ntf.end_pipeline.assert_called_with(1)
        mocked_finish.assert_called_with(c.DATASET_PROCESSED_FAIL)

    @patched_finish
    @patch(ppath('MostRecentProc.retrieve_entity'))
    def test_abort(self, mocked_retrieve_entity, mocked_finish):
        self.dataset.abort()
        mocked_finish.assert_called_with(c.DATASET_ABORTED)

    @patch(ppath('MostRecentProc.change_status'))
    @patch(ppath('MostRecentProc.retrieve_entity'))
    def test_reset(self, mocked_retrieve_entity, mocked_change_status):
        self.dataset.reset()
        mocked_change_status.assert_called_with(c.DATASET_REPROCESS)

    @patch(ppath('MostRecentProc.start_stage'))
    def test_start_stage(self, mocked_start_stage):
        self.dataset.start_stage('a_stage')
        mocked_start_stage.assert_called_with('a_stage')
        self.dataset.ntf.start_stage.assert_called_with('a_stage')

    @patch(ppath('MostRecentProc.end_stage'))
    def test_end_stage(self, mocked_end_stage):
        self.dataset.end_stage('a_stage')
        mocked_end_stage.assert_called_with('a_stage', 0)
        self.dataset.ntf.end_stage.assert_called_with('a_stage', 0)

    def test_str(self):
        with patched_stages:
            assert str(self.dataset) == 'test_dataset -- this, that, other'

    def setup_dataset(self):
        with patched_get():
            self.dataset = _TestDataset(
                'test_dataset',
                {'date_started': 'now', 'dataset_name': 'None', 'dataset_type': 'None'}
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
                d = RunDataset(d_name, os.path.join(self.base_dir, d_name))
                assert d._is_ready() == rta_complete

    def test_dataset_status(self):
        super().test_dataset_status()
        del self.dataset.most_recent_proc.entity['status']
        assert not self.dataset._is_ready()
        assert self.dataset.dataset_status == c.DATASET_NEW
        os.mkdir(os.path.join(self.base_dir, self.dataset.name))
        touch(os.path.join(self.base_dir, self.dataset.name, 'RTAComplete.txt'))
        assert self.dataset._is_ready()
        assert self.dataset.dataset_status == c.DATASET_READY

    def setup_dataset(self):
        self.dataset = RunDataset(
            'test_dataset',
            os.path.join(self.base_dir, 'test_dataset'),
            most_recent_proc={'date_started': 'now', 'dataset_name': 'None', 'dataset_type': 'None'}
        )


class TestSampleDataset(TestDataset):
    def test_dataset_status(self):
        with patched_expected_yield():
            super().test_dataset_status()

    @patch(ppath('MostRecentProc.change_status'))
    def test_force(self, mocked_change_status):
        self.dataset.force()
        mocked_change_status.assert_called_with(c.DATASET_FORCE_READY)

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
        expected_str = 'test_dataset -- this, that, other  (480 / 1000000000  from a_run_id, another_run_id) (non useable in a_run_id, another_run_id)'
        self.dataset._data_threshold = None
        with patched_get(self.dataset.run_elements), patched_expected_yield(), patched_stages:
            print(expected_str)
            print(str(self.dataset))
            assert str(self.dataset) == expected_str

    def setup_dataset(self):
        with patched_get():
            self.dataset = SampleDataset(
                'test_dataset',
                most_recent_proc={'date_started': 'now', 'dataset_name': 'None', 'dataset_type': 'None'}
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


class TestMostRecentProc(TestAnalysisDriver):
    def setUp(self):
        with patched_get():
            self.proc = MostRecentProc('test', 'test')
        self.proc._entity = fake_proc.copy()

    def test_rest_entity_pre_existing(self):
        assert self.proc.entity == fake_proc

    @patched_get()
    @patched_initialise
    def test_rest_entity_not_pre_existing(self, mocked_initialise, mocked_get):
        self.proc._entity = None
        with patched_datetime():
            x = self.proc.entity
            assert x == fake_proc
            mocked_get.assert_called_with(
                'analysis_driver_procs',
                where={'dataset_type': 'test', 'dataset_name': 'test'},
                sort='-_created'
            )

        self.proc._entity = None
        mocked_get.return_value = []
        y = self.proc.entity
        assert y == {}
        mocked_get.assert_called_with(
            'analysis_driver_procs',
            where={'dataset_type': 'test', 'dataset_name': 'test'},
            sort='-_created'
        )
        assert mocked_initialise.call_count == 0

    @patched_post
    @patched_patch
    def test_initialise_entity(self, mocked_patch, mocked_post):
        self.proc._entity = None
        with patched_datetime():
            self.proc.initialise_entity()
        mocked_post.assert_called_with(
            'analysis_driver_procs',
            {'proc_id': 'test_test_now', 'dataset_type': 'test', 'dataset_name': 'test'}
        )
        mocked_patch.assert_called_with(
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
        self.proc.sync()
        assert self.proc._entity == {
            'proc_id': 'test_test_now',
            'dataset_type': 'test',
            'dataset_name': 'test'
        }

        self.proc.entity.update({'this': 'that'})
        self.proc.sync()
        mocked_patch.assert_called_with(
            'analysis_driver_procs',
            {
                'this': 'that',
                'dataset_type': 'test',
                'dataset_name': 'test'
            },
            id_field='proc_id',
            element_id='test_test_now'
        )
        assert self.proc._entity == {
            'proc_id': 'test_test_now',
            'dataset_type': 'test',
            'dataset_name': 'test',
            'this': 'that'
        }

    @patch(ppath('MostRecentProc.sync'))
    def test_update_entity(self, mocked_sync):
        self.proc.update_entity(other='another')
        assert self.proc.entity == {
            'proc_id': 'test_test_now',
            'dataset_type': 'test',
            'dataset_name': 'test',
            'other': 'another'
        }
        mocked_sync.assert_called()

    @patched_update
    def test_start(self, mocked_update):
        with patched_pid:
            self.proc.start()
        mocked_update.assert_called_with(status=c.DATASET_PROCESSING, pid=1)

    @patched_update
    def test_finish(self, mocked_update):
        with patched_datetime():
            self.proc.finish(c.DATASET_PROCESSED_SUCCESS)
        mocked_update.assert_called_with(status=c.DATASET_PROCESSED_SUCCESS, pid=None, end_date='now')

    @patched_update
    @patched_post
    def test_start_stage(self, mocked_post, mocked_update):
        with patched_datetime('then'), patched_stage_id('a_stage_id'):
            self.proc.start_stage('test_stage')

        mocked_update.assert_called_with(stages=['a_stage_id'])
        mocked_post.assert_called_with(
            'analysis_driver_stages',
            {'date_started': 'then', 'stage_name': 'test_stage', 'analysis_driver_proc': 'test_test_now'}
        )

        self.proc.entity['stages'] = ['a_stage_id']
        with patched_datetime('later'), patched_stage_id('another_stage_id'):
            self.proc.start_stage('another_stage')

        mocked_update.assert_called_with(stages=['a_stage_id', 'another_stage_id'])
        mocked_post.assert_called_with(
            'analysis_driver_stages',
            {'date_started': 'later', 'stage_name': 'another_stage', 'analysis_driver_proc': 'test_test_now'}
        )

    @patched_patch
    @patched_datetime()
    @patched_stage_id('a_stage_id')
    def test_end_stage(self, mocked_stage_id, mocked_now, mocked_patch):
        self.proc.end_stage('this')

        mocked_stage_id.assert_called_with('this')
        mocked_patch.assert_called_with(
            'analysis_driver_stages', {'date_finished': 'now', 'exit_status': 0}, '_id', 'a_stage_id'
        )

    @patch(ppath('rest_communication', 'get_document'), return_value={'_id': 'a_stage_id'})
    def test_stage_id(self, mocked_get):
        assert self.proc._stage_id('a_stage_name')
        mocked_get.assert_called_with(
            'analysis_driver_stages',
            where={'analysis_driver_proc': 'test_test_now', 'stage_name': 'a_stage_name'}
        )
