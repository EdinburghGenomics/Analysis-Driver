import os
from os.path import join as p_join
import shutil
from tests.test_analysisdriver import TestAnalysisDriver
from unittest.mock import patch, Mock
from bin import delete_data
from analysis_driver.util import find_files, find_fastqs


class FakeExecutor(Mock):
    @staticmethod
    def join():
        pass


class TestDeleter(TestAnalysisDriver):
    @property
    def assets_deletion(self):
        return p_join(self.assets_path, 'data_deletion')

    @staticmethod
    def touch(file_path):
        open(file_path, 'w').close()

    def setUp(self):
        self.deleter = delete_data.Deleter(self.assets_deletion)

    @patch('bin.delete_data.executor.execute', return_value=FakeExecutor())
    def test_execute(self, mocked_execute):
        self.deleter._execute('a test command')
        mocked_execute.assert_called_with(['a test command'], 'local')


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


class TestRawDataDeleter(TestDeleter):
    def _setup_run(self, run_id, deletable_sub_dirs):
        for d in deletable_sub_dirs + ('Stats', 'InterOp', 'RTAComplete.txt'):
            os.makedirs(p_join(self.assets_deletion, 'raw', run_id, d), exist_ok=True)

    def setUp(self):
        self.original_input_dir = delete_data.cfg['input_dir']
        delete_data.cfg.content['input_dir'] = p_join(self.assets_deletion, 'raw')
        self.deleter = delete_data.RawDataDeleter(self.assets_deletion)
        self._setup_run('deletable_run', self.deleter.deletable_sub_dirs)
        os.makedirs(p_join(self.assets_deletion, 'archive'), exist_ok=True)

    def tearDown(self):
        delete_data.cfg.content['input_dir'] = self.original_input_dir
        shutil.rmtree(p_join(self.assets_deletion, 'raw', 'deletable_run'), ignore_errors=True)
        for d in os.listdir(p_join(self.assets_deletion, 'archive')):
            shutil.rmtree(p_join(self.assets_deletion, 'archive', d))

        deletion_script = p_join(self.assets_deletion, 'data_deletion.pbs')
        if os.path.isfile(deletion_script):
            os.remove(deletion_script)

    def test_deletable_runs(self):
        patch_target = 'analysis_driver.rest_communication.get_documents'
        expected_rest_query = {
            'max_results': 100,
            'embedded': {'run_elements': 1, 'analysis_driver_procs': 1},
            'aggregate': True,
            'sort': 'run_id'
        }

        with patch(patch_target, return_value=fake_run_elements_no_procs) as p:
            runs = self.deleter.deletable_runs()
            p.assert_called_with('runs', depaginate=True, **expected_rest_query)
            assert runs == []

        with patch(patch_target, return_value=fake_run_elements_procs_running) as p:
            runs = self.deleter.deletable_runs()
            p.assert_called_with('runs', depaginate=True, **expected_rest_query)
            assert runs == []

        with patch(patch_target, return_value=fake_run_elements_procs_complete) as p:
            runs = self.deleter.deletable_runs()
            p.assert_called_with('runs', depaginate=True, **expected_rest_query)
            assert runs == fake_run_elements_procs_complete[1:]

    def test_setup_run_for_deletion(self):
        deletion_dir = p_join(self.assets_deletion, 'raw', '.deletion')
        subdirs_to_delete = self.deleter._setup_run_for_deletion(
            'deletable_run',
            deletion_dir=deletion_dir
        )
        self.compare_lists(subdirs_to_delete, self.deleter.deletable_sub_dirs)
        shutil.rmtree(deletion_dir)

    def test_setup_runs_for_deletion(self):
        with patched_deletable_runs:
            deletion_dir = self.deleter.setup_runs_for_deletion(self.deleter.deletable_runs())
            assert os.listdir(deletion_dir) == ['deletable_run']
        shutil.rmtree(deletion_dir)

    def test_delete_runs(self):
        with patched_deletable_runs:
            deletion_dir = self.deleter.setup_runs_for_deletion(self.deleter.deletable_runs())
            self.deleter.delete_dir(deletion_dir)
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
            'proc_id',
            'most_recent_proc'
        )

    def test_archive_run(self):
        run_id = 'run_to_archive'
        raw_dir = p_join(self.assets_deletion, 'raw', run_id)
        self._setup_run(run_id, deletable_sub_dirs=())

        self.deleter.archive_run(run_id)
        assert not os.path.isdir(raw_dir)
        archived_run = p_join(self.assets_deletion, 'archive', run_id)
        assert os.path.isdir(archived_run)
        self.compare_lists(
            os.listdir(archived_run),
            ['InterOp', 'RTAComplete.txt', 'Stats']
        )
        shutil.rmtree(archived_run)

    @patched_patch_entry
    @patched_deletable_runs
    def test_run_deletion(self, mocked_deletable_runs, mocked_patch):
        self._setup_run('non_deletable_run', self.deleter.deletable_sub_dirs)
        del_dir = p_join(self.assets_deletion, 'raw', 'deletable_run')
        non_del_dir = p_join(self.assets_deletion, 'raw', 'non_deletable_run')
        self.compare_lists(
            os.listdir(p_join(self.assets_deletion, 'raw')),
            ['deletable_run', 'non_deletable_run']
        )
        for d in (del_dir, non_del_dir):
            self.compare_lists(
                os.listdir(d),
                list(self.deleter.deletable_sub_dirs) + ['Stats', 'InterOp', 'RTAComplete.txt']
            )
        self.deleter.delete_data()
        assert os.path.isdir(non_del_dir)
        self.compare_lists(
            os.listdir(non_del_dir),
            list(self.deleter.deletable_sub_dirs) + ['Stats', 'InterOp', 'RTAComplete.txt']
        )
        self.compare_lists(
            os.listdir(p_join(self.assets_deletion, 'archive', 'deletable_run')),
            ['Stats', 'InterOp', 'RTAComplete.txt']
        )
        mocked_patch.assert_called_with(
            'analysis_driver_procs',
            {'status': 'deleted'},
            'proc_id',
            'most_recent_proc'
        )
        assert mocked_deletable_runs.call_count == 1


patched_clarity_get_samples = patch(
    'bin.delete_data.clarity.get_released_samples', return_value=['deletable_sample', 'deletable_sample_2']
)

patched_deletable_samples = patch(
    'bin.delete_data.rest_communication.get_documents',
    return_value=[{'sample_id': 'deletable_sample'}]
)


class TestFastqDeleter(TestDeleter):
    @property
    def fake_deletion_record(self):
        return delete_data._FastqDeletionRecord(
            run_element={
                'run_id': 'a_run',
                'project_id': 'a_project',
                'sample_id': 'deletable_sample',
                'lane': '1'
            },
            fastqs=find_fastqs(
                p_join(self.assets_deletion, 'fastqs', 'a_run', 'fastq'),
                'a_project',
                'deletable_sample',
                lane=1
            )
        )

    def _setup_fastqs(self, run_id, project_id, sample_id):
        fastq_dir = p_join(self.assets_deletion, 'fastqs', run_id, 'fastq', project_id, sample_id)
        os.makedirs(fastq_dir, exist_ok=True)
        for lane in range(8):
            for read in ('1', '2'):
                for file_ext in ('fastq.gz', 'fastqc.html'):
                    self.touch(p_join(fastq_dir, 'fastq_L00%s_R%s.%s' % (str(lane + 1), read, file_ext)))

    def setUp(self):
        self.deleter = delete_data.FastqDeleter(self.assets_deletion)
        for run_id in ('a_run', 'another_run'):
            for sample_id in ('deletable_sample', 'non_deletable_sample'):
                self._setup_fastqs(run_id, 'a_project', sample_id)

    def tearDown(self):
        for r in ('a_run', 'another_run'):
            shutil.rmtree(p_join(self.assets_deletion, 'fastqs', r), ignore_errors=True)

        deletion_script = p_join(self.assets_deletion, 'data_deletion.pbs')
        if os.path.isfile(deletion_script):
            os.remove(deletion_script)

        for tmpdir in find_files(self.assets_deletion, 'fastqs', '.data_deletion_*'):
            shutil.rmtree(tmpdir)

    def test_samples_released_in_lims(self):
        with patched_clarity_get_samples:
            assert self.deleter.samples_released_in_lims == {'deletable_sample', 'deletable_sample_2'}

    def test_samples_released_in_app(self):
        with patched_deletable_samples:
            self.compare_lists(self.deleter.samples_released_in_app, ['deletable_sample'])

    def test_find_fastqs_for_run_element(self):
        run_element = {
            'run_id': 'a_run',
            'project_id': 'a_project',
            'sample_id': 'deletable_sample',
            'lane': '1'
        }
        fqs = self.deleter.find_fastqs_for_run_element(run_element)
        assert len(fqs) == 2

        obs = [os.path.basename(f) for f in fqs]
        expected_fqs = find_fastqs(
            p_join(self.assets_deletion, 'fastqs', 'a_run', 'fastq'),
            'a_project',
            'deletable_sample',
            lane=1
        )
        exp = [os.path.basename(f) for f in expected_fqs]
        self.compare_lists(obs, exp)

    def test_setup_record_for_deletion(self):
        deletion_dir = p_join(self.assets_deletion, 'fastqs', '.test_setup_record')
        os.makedirs(deletion_dir, exist_ok=True)
        e = {
            'run_id': 'a_run',
            'project_id': 'a_project',
            'sample_id': 'deletable_sample',
            'lane': '1'
        }
        fqs = self.deleter.find_fastqs_for_run_element(e)
        record = delete_data._FastqDeletionRecord(e, fqs)
        self.deleter._setup_record_for_deletion(deletion_dir, record)

        self.compare_lists(os.listdir(deletion_dir), ['a_run'])
        self.compare_lists(os.listdir(p_join(deletion_dir, 'a_run', 'fastq', 'a_project')), ['deletable_sample'])
        self.compare_lists(
            os.listdir(p_join(deletion_dir, 'a_run', 'fastq', 'a_project', 'deletable_sample')),
            [os.path.basename(fq) for fq in fqs]
        )
        shutil.rmtree(deletion_dir)

    def setup_deletion_records(self):  # test set() intersection between lims and app
        with patched_clarity_get_samples, patched_deletable_samples:
            records = self.deleter.setup_deletion_records()
            assert len(records) == 1 and records[0].sample_id == 'deletable_sample'

    def test_setup_fastqs_for_deletion(self):
        records = [self.fake_deletion_record]
        with patched_clarity_get_samples:
            self.deleter.setup_fastqs_for_deletion(records)

    def test_delete_data(self):
        run_elements = []
        for run in ('a_run', 'another_run'):
            for lane in range(8):
                e = {
                    'run_id': run,
                    'project_id': 'a_project',
                    'sample_id': 'deletable_sample',
                    'lane': str(lane + 1)
                }
                run_elements.append(e)

        patched_app = patch(
            'bin.delete_data.FastqDeleter.samples_released_in_app',
            new={'deletable_sample'}
        )
        patched_run_elements = patch(
            'bin.delete_data.rest_communication.get_documents',
            return_value=run_elements
        )
        patched_mark_sample = patch(
            'bin.delete_data.FastqDeleter.mark_sample_fastqs_as_deleted',
        )
        basenames = []
        for lane in range(8):
            for read in ('1', '2'):
                for file_ext in ('fastq.gz', 'fastqc.html'):
                    basenames.append('fastq_L00%s_R%s.%s' % (lane + 1, read, file_ext))

        with patched_clarity_get_samples, patched_app, patched_run_elements, patched_mark_sample:
            self.compare_lists(
                [os.path.basename(f) for f in find_files(self.assets_deletion, 'fastqs', '*run*')],
                ['a_run', 'another_run']
            )
            for run in ('a_run', 'another_run'):
                self.compare_lists(
                    os.listdir(p_join(self.assets_deletion, 'fastqs', run, 'fastq', 'a_project', 'deletable_sample')),
                    basenames
                )
            self.deleter.delete_data()

            for run in ('a_run', 'another_run'):
                self.compare_lists(
                    os.listdir(p_join(self.assets_deletion, 'fastqs', run, 'fastq', 'a_project', 'deletable_sample')),
                    [f for f in basenames if f.endswith('fastqc.html')]
                )
