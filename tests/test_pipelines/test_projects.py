import os
from shutil import rmtree
from unittest.mock import patch
from tests.test_analysisdriver import TestAnalysisDriver, NamedMock
from analysis_driver.pipelines import projects

test_projects = os.path.join(TestAnalysisDriver.assets_path, 'test_projects')
relatedness_outfiles = os.path.join(test_projects, 'relatedness_outfiles')


class TestProjects(TestAnalysisDriver):
    def test_build_pipeline(self):
        project_id = 'test_dataset'
        two_sample_dataset = NamedMock(
            real_name=project_id,
            species='Homo sapiens',
            samples_processed=[
                {
                    'sample_id': '10015AT0004', 'user_sample_id': 'test_user_sample1',
                    'aggregated': {'most_recent_proc': {'pipeline_used': {'name': 'bcbio'}}}
                },
                {
                    'sample_id': '10015AT0003', 'user_sample_id': 'test_user_sample2',
                    'aggregated': {'most_recent_proc': {'pipeline_used': {'name': 'bcbio'}}}
                }
            ]
        )
        pipeline = projects.build_pipeline(two_sample_dataset)
        # get a real pipeline
        assert len(pipeline.previous_stages) > 0
        one_sample_dataset = NamedMock(
            real_name=project_id,
            species='Homo sapiens',
            samples_processed=[{'sample_id': '10015AT0004', 'user_sample_id': 'test_user_sample1'}]
        )
        pipeline = projects.build_pipeline(one_sample_dataset)
        # get a single step pipeline
        assert len(pipeline.previous_stages) == 0

        multi_non_human_samples_dataset = NamedMock(
            real_name=project_id,
            species='Thingymy',
            samples_processed=[
                {
                    'sample_id': 'sample1', 'user_sample_id': 'user_sample1',
                    'aggregated': {'most_recent_proc': {'pipeline_used': {'name': 'qc'}}}
                },
                {
                    'sample_id': 'sample2', 'user_sample_id': 'user_sample2',
                    'aggregated': {'most_recent_proc': {'pipeline_used': {'name': 'qc'}}}
                },
                {
                    'sample_id': 'sample3', 'user_sample_id': 'user_sample3',
                    'aggregated': {'most_recent_proc': {'pipeline_used': {'name': 'qc'}}}
                },
            ]
        )
        pipeline = projects.build_pipeline(multi_non_human_samples_dataset)
        # get a single step pipeline
        assert len(pipeline.previous_stages) == 0


class TestMD5Sum(TestAnalysisDriver):
    def setUp(self):
        os.makedirs(relatedness_outfiles, exist_ok=True)
        self.filename = os.path.join(relatedness_outfiles, 'an_output_file.txt')
        open(self.filename, 'w').close()

    def tearDown(self):
        os.remove(self.filename)
        rmtree(relatedness_outfiles)

    @patch.object(projects.MD5Sum, 'job_dir', new=test_projects)
    @patch('egcg_core.executor.execute')
    def test_run(self, mocked_execute):
        dataset = NamedMock(
            real_name='test_dataset',
            samples_processed=[
                {'sample_id': '10015AT0004', 'user_sample_id': 'test_user_sample1'},
                {'sample_id': '10015AT0003', 'user_sample_id': 'test_user_sample2'}
            ]
        )
        md5 = projects.MD5Sum(dataset=dataset)
        with patch('os.makedirs'):
            md5._run()

        mocked_execute.assert_called_with(
            'path/to/md5sum %s > %s' % (self.filename, self.filename + '.md5'),
            cpus=1,
            mem=2,
            job_name='md5sum',
            log_commands=False,
            working_dir=test_projects
        )


class TestOutput(TestAnalysisDriver):
    def setUp(self):
        dataset = NamedMock(
            real_name='test_dataset',
            samples_processed=[
                {'sample_id': '10015AT0004', 'user_sample_id': 'test_user_sample1'},
                {'sample_id': '10015AT0003', 'user_sample_id': 'test_user_sample2'}
            ]
        )
        self.o = projects.Output(dataset=dataset)

    @patch.object(projects.toolset, 'write_to_yaml')
    @patch.object(projects, 'create_output_links')
    @patch.object(projects, 'output_data_and_archive')
    @patch.object(projects, 'output_file_config')
    def test_run(self, mocked_outfile_config, mocked_output_archive, mocked_output_links, mocked_write):
        with patch('analysis_driver.segmentation.BasicStage.job_dir', new=test_projects):
            self.o._run()
            mocked_output_archive.assert_called_with(relatedness_outfiles, '/path/to/input/dir/test_dataset')
            mocked_outfile_config.set_pipeline_type.assert_called_with('project_process')
            mocked_write.assert_called_with(os.path.join(relatedness_outfiles, 'program_versions.yaml'))
            mocked_output_links.assert_called_with(
                test_projects,
                mocked_outfile_config,
                relatedness_outfiles,
                project_id='test_dataset'
            )
