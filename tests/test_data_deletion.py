__author__ = 'mwham'
import os
from os.path import join
import shutil
from tests.test_analysisdriver import TestAnalysisDriver
from unittest.mock import patch, Mock
from bin import delete_data


finished_proc = {'status': 'finished'}
aborted_proc = {'status': 'aborted'}
running_proc = {'status': 'running'}


# none of these should be deletable
fake_run_elements_no_procs = [
    {'run_element_id': 'unreviewed_run_element', 'run_id': 'unreviewed_run', 'review_statuses': ['not reviewed']},
    {'run_element_id': 'passed_run_element', 'run_id': 'passed_run', 'review_statuses': ['pass']},
    {'run_element_id': 'failed_run_element', 'run_id': 'failed_run', 'review_statuses': ['fail']}
]


fake_run_elements_procs_running = [
    {
        # not deletable
        'run_element_id': 'unreviewed_run_element',
        'run_id': 'unreviewed_run',
        'review_statuses': ['not reviewed'],
        'analysis_driver_procs': [running_proc]
    },
    {
        # not deletable
        'run_element_id': 'partially_unreviewed_run_element',
        'run_id': 'partially_unreviewed_run',
        'review_statuses': ['not reviewed', 'pass'],
        'analysis_driver_procs': [running_proc]
    },
    {
        # not deletable
        'run_element_id': 'passed_run_element',
        'run_id': 'passed_run',
        'review_statuses': ['pass'],
        'analysis_driver_procs': [running_proc]
    },
    {
        # not deletable
        'run_element_id': 'failed_run_element',
        'run_id': 'failed_run',
        'review_statuses': ['fail'],
        'analysis_driver_procs': [running_proc]
    }
]


fake_run_elements_procs_complete = [
    {
        # not deletable
        'run_element_id': 'unreviewed_run_element',
        'run_id': 'unreviewed_run',
        'review_statuses': ['not reviewed', 'pass'],
        'analysis_driver_procs': [running_proc, finished_proc]
    },
    {
        # deletable
        'run_element_id': 'passed_run_element',
        'run_id': 'passed_run',
        'review_statuses': ['pass'],
        'analysis_driver_procs': [running_proc, finished_proc]
    },
    {
        # deletable
        'run_element_id': 'failed_run_element',
        'run_id': 'failed_run',
        'review_statuses': ['fail'],
        'analysis_driver_procs': [running_proc, aborted_proc]
    }

]


class FakeExecutor(Mock):
    @staticmethod
    def join():
        pass


patched_executor = patch('bin.delete_data.executor.execute', return_value=FakeExecutor())
patched_deletable_runs = patch(
    'bin.delete_data.RawDataDeleter.deletable_runs',
    return_value=[{'run_id': 'deletable_run'}]
)


class TestDeleter(TestAnalysisDriver):
    @property
    def assets_deletion(self):
        return join(self.assets_path, 'data_deletion')

    def setUp(self):
        self.deleter = delete_data.Deleter(self.assets_deletion)

    @patched_executor
    def test_execute(self, mocked_execute):
        self.deleter._execute('a test command')
        mocked_execute.assert_called_with(['a test command'], 'local')


class TestRawDataDeleter(TestDeleter):
    expected_rest_query = (
        'runs', 'embedded={"run_elements":1,"analysis_driver_procs":1}', 'aggregate=True', 'max_results=1000'
    )

    def _setup_run(self, run_id, deletable_sub_dirs):
        for d in deletable_sub_dirs + ('Stats', 'InterOp', 'RTAComplete.txt'):
            os.makedirs(join(self.assets_deletion, 'raw', run_id, d))

    def setUp(self):
        self.original_input_dir = delete_data.cfg['input_dir']
        delete_data.cfg.content['input_dir'] = join(self.assets_deletion, 'raw')
        self.deleter = delete_data.RawDataDeleter(self.assets_deletion)
        self._setup_run('deletable_run', self.deleter.deletable_sub_dirs)

    def tearDown(self):
        delete_data.cfg.content['input_dir'] = self.original_input_dir
        shutil.rmtree(join(self.assets_deletion, 'raw', 'deletable_run'))

    def test_deletable_runs(self):
        patch_target = 'bin.delete_data.query_api'

        with patch(patch_target, return_value=fake_run_elements_no_procs) as p:
            runs = self.deleter.deletable_runs()
            p.assert_called_with(*self.expected_rest_query)
            assert runs == []

        with patch(patch_target, return_value=fake_run_elements_procs_running) as p:
            runs = self.deleter.deletable_runs()
            p.assert_called_with(*self.expected_rest_query)
            assert runs == []

        with patch(patch_target, return_value=fake_run_elements_procs_complete) as p:
            runs = self.deleter.deletable_runs()
            p.assert_called_with(*self.expected_rest_query)
            assert runs == fake_run_elements_procs_complete[1:]

    def test_setup_run_for_deletion(self):
        deletion_dir = join(self.assets_deletion, 'raw', '.deletion')
        subdirs_to_delete = self.deleter._setup_run_for_deletion(
            'deletable_run',
            deletion_dir=deletion_dir
        )
        self.compare_lists(subdirs_to_delete, self.deleter.deletable_sub_dirs)
        shutil.rmtree(deletion_dir)

    def test_setup_runs_for_deletion(self):
        with patched_deletable_runs:
            deletion_dir = self.deleter.setup_runs_for_deletion()
            assert os.listdir(deletion_dir) == ['deletable_run']
        shutil.rmtree(deletion_dir)

    def test_delete_runs(self):
        with patched_deletable_runs:
            self.deleter.delete_runs(self.deleter.setup_runs_for_deletion())
