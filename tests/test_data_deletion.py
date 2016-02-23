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


patched_patch_entry = patch('bin.delete_data.rest_communication.patch_entry')
patched_deletable_runs = patch(
    'bin.delete_data.RawDataDeleter.deletable_runs',
    return_value=[
        {
            'run_id': 'deletable_run',
            'analysis_driver_procs': [
                {'proc_id': 'first_proc'},
                {'proc_id': 'most_recent_proc'}
            ]
        }
    ]
)


class TestDeleter(TestAnalysisDriver):
    @property
    def assets_deletion(self):
        return join(self.assets_path, 'data_deletion')

    def setUp(self):
        self.deleter = delete_data.Deleter(self.assets_deletion)

    @patch('bin.delete_data.executor.execute', return_value=FakeExecutor())
    def test_execute(self, mocked_execute):
        self.deleter._execute('a test command')
        mocked_execute.assert_called_with(['a test command'], 'local')


class TestRawDataDeleter(TestDeleter):
    expected_rest_query = (
        'runs', 'embedded={"run_elements":1,"analysis_driver_procs":1}', 'aggregate=True', 'max_results=50'
    )

    def _setup_run(self, run_id, deletable_sub_dirs):
        for d in deletable_sub_dirs + ('Stats', 'InterOp', 'RTAComplete.txt'):
            os.makedirs(join(self.assets_deletion, 'raw', run_id, d), exist_ok=True)

    def setUp(self):
        self.original_input_dir = delete_data.cfg['input_dir']
        delete_data.cfg.content['input_dir'] = join(self.assets_deletion, 'raw')
        self.deleter = delete_data.RawDataDeleter(self.assets_deletion)
        self._setup_run('deletable_run', self.deleter.deletable_sub_dirs)
        os.makedirs(join(self.assets_deletion, 'archive'), exist_ok=True)

    def tearDown(self):
        delete_data.cfg.content['input_dir'] = self.original_input_dir
        shutil.rmtree(join(self.assets_deletion, 'raw', 'deletable_run'), ignore_errors=True)
        for d in os.listdir(join(self.assets_deletion, 'archive')):
            shutil.rmtree(join(join(self.assets_deletion, 'archive', d)))

        deletion_script = join(self.assets_deletion, 'data_deletion.pbs')
        if os.path.isfile(deletion_script):
            os.remove(deletion_script)

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
            deletion_dir = self.deleter.setup_runs_for_deletion()
            self.deleter.delete_runs(deletion_dir)
            assert not os.path.isdir(deletion_dir)

    @patched_patch_entry
    def test_mark_run_as_deleted(self, mocked_patch):
        run_object = {
            'run_id': 'a_run',
            'analysis_driver_procs': [
                {'proc_id': 'first_proc'},
                {'proc_id': 'most_recent_proc'}
            ]
        }
        self.deleter.mark_run_as_deleted(run_object)
        mocked_patch.assert_called_with(
            'analysis_driver_procs',
            {'status': 'deleted'},
            proc_id='most_recent_proc'
        )

    def test_archive_run(self):
        run_id = 'run_to_archive'
        raw_dir = join(self.assets_deletion, 'raw', run_id)
        self._setup_run(run_id, deletable_sub_dirs=())

        self.deleter.archive_run(run_id)
        assert not os.path.isdir(raw_dir)
        archived_run = join(self.assets_deletion, 'archive', run_id)
        assert os.path.isdir(archived_run)
        assert os.listdir(archived_run) == ['InterOp', 'RTAComplete.txt', 'Stats']
        shutil.rmtree(archived_run)

    @patched_patch_entry
    @patched_deletable_runs
    def test_run_deletion(self, mocked_deletable_runs, mocked_patch):
        self._setup_run('non_deletable_run', self.deleter.deletable_sub_dirs)
        del_dir = join(self.assets_deletion, 'raw', 'deletable_run')
        non_del_dir = join(self.assets_deletion, 'raw', 'non_deletable_run')
        self.compare_lists(
            os.listdir(join(self.assets_deletion, 'raw')),
            ['deletable_run', 'non_deletable_run']
        )
        for d in (del_dir, non_del_dir):
            self.compare_lists(
                os.listdir(d),
                list(self.deleter.deletable_sub_dirs) + ['Stats', 'InterOp', 'RTAComplete.txt']
            )
        self.deleter.run_deletion()
        assert os.path.isdir(non_del_dir)
        self.compare_lists(
            os.listdir(non_del_dir),
            list(self.deleter.deletable_sub_dirs) + ['Stats', 'InterOp', 'RTAComplete.txt']
        )
        self.compare_lists(
            os.listdir(join(self.assets_deletion, 'archive', 'deletable_run')),
            ['Stats', 'InterOp', 'RTAComplete.txt']
        )
        mocked_patch.assert_called_with(
            'analysis_driver_procs',
            {'status': 'deleted'},
            proc_id='most_recent_proc'
        )
        assert mocked_deletable_runs.call_count == 2
