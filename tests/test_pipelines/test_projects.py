import os
from unittest.mock import patch
from tests.test_analysisdriver import TestAnalysisDriver, NamedMock
from analysis_driver.pipelines import projects
from analysis_driver.exceptions import PipelineError


class TestProjects(TestAnalysisDriver):
    def setUp(self):
        self.project_id = 'test_dataset'
        self.dataset = NamedMock(real_name=self.project_id,
                                 samples_processed=[{'sample_id': '10015AT0004', 'user_sample_id': 'test_user_sample1'},
                                                    {'sample_id': '10015AT0003', 'user_sample_id': 'test_user_sample2'}],
                                 name='10015AT0004',
                                 species='Homo sapiens')

    def test_build_pipeline(self):
        with patch('analysis_driver.pipelines.projects.delivery_source', return_value=os.path.join(self.assets_path, 'test_projects')):
            projects.build_pipeline(self.dataset)

        with patch('analysis_driver.pipelines.projects.delivery_source', return_value=os.path.join(self.assets_path, 'nonexistent_project')):
            with self.assertRaises(PipelineError):
                projects.build_pipeline(self.dataset)

        self.dataset = NamedMock(real_name=self.project_id,
                                 samples_processed=[{'sample_id': '10015AT0004', 'user_sample_id': 'test_user_sample1'}],
                                 name='10015AT0004',
                                 species='Homo sapiens')

        with patch('analysis_driver.pipelines.projects.delivery_source', return_value=os.path.join(self.assets_path, 'test_projects')):
            with self.assertRaises(PipelineError):
                projects.build_pipeline(self.dataset)


class TestMD5Sum(TestAnalysisDriver):
    def setUp(self):
        self.project_id = 'test_dataset'
        dataset = NamedMock(real_name=self.project_id, samples_processed=[{'sample_id': '10015AT0004', 'user_sample_id': 'test_user_sample1'},
                                                                               {'sample_id': '10015AT0003', 'user_sample_id': 'test_user_sample2'}])
        self.md5 = projects.MD5Sum(dataset=dataset)


    @patch('egcg_core.executor.execute')
    def test_run(self, mocked_execute):
        with patch('analysis_driver.segmentation.BasicStage.job_dir', new=os.path.join(self.assets_path, 'test_projects')):
            self.md5._run()


        mocked_execute.assert_called_once_with(cpus=1,
                                               job_name='md5sum',
                                               log_commands=False,
                                               mem=2,
                                               working_dir='/Users/katie/PycharmProjects/Analysis-Driver/tests/assets/test_projects')


class TestOutput(TestAnalysisDriver):
    def setUp(self):
        self.project_id = 'test_dataset'
        dataset = NamedMock(real_name=self.project_id, samples_processed=[{'sample_id': '10015AT0004', 'user_sample_id': 'test_user_sample1'},
                                                                               {'sample_id': '10015AT0003', 'user_sample_id': 'test_user_sample2'}])
        self.o = projects.Output(dataset=dataset)

    @patch('analysis_driver.pipelines.projects.create_output_links')
    @patch('analysis_driver.pipelines.projects.output_data_and_archive')
    @patch('analysis_driver.pipelines.projects.OutputFileConfiguration', return_value='OutfileConfig')
    def test_run(self, mocked_outfile_config, mocked_output_archive, mocked_output_links):
        with patch('analysis_driver.segmentation.BasicStage.job_dir', new=os.path.join(self.assets_path, 'test_projects')):
            self.o._run()
            assert mocked_output_archive.called_with('/Users/katie/PycharmProjects/Analysis-Driver/tests/assets/test_projects/relatedness_outfiles',
                                                     '/path/to/input/dir/test_dataset')
            assert mocked_output_links.called_with('/Users/katie/PycharmProjects/Analysis-Driver/tests/assets/test_projects',
                                                   'OutfileConfig',
                                                   '/Users/katie/PycharmProjects/Analysis-Driver/tests/assets/test_projects/relatedness_outfiles')
