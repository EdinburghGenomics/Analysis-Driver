from os.path import join
from tests.test_analysisdriver import TestAnalysisDriver, NamedMock


class QCTester(TestAnalysisDriver):
    def setUp(self):
        self.sample_id = 'test_sample'
        self.run_id = 'test_run'
        self.dataset = NamedMock(real_name=self.sample_id, data_threshold=None)
        self.run_dataset = NamedMock(real_name=self.run_id)

    @staticmethod
    def fake_find_file(*path_parts):
        return join(*path_parts)
