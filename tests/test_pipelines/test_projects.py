import os
from unittest.mock import Mock, patch, PropertyMock
from tests.test_analysisdriver import TestAnalysisDriver, NamedMock
from analysis_driver.pipelines import projects


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



#class TestMD5Sum(TestProjects):
#    pass

#class Output(TestProjects):
#    pass
