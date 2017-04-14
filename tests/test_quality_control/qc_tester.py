from tests.test_analysisdriver import TestAnalysisDriver, NamedMock


class QCTester(TestAnalysisDriver):
    def setUp(self):
        self.dataset = NamedMock(real_name='test_sample', data_threshold=None)
        self.run_dataset = NamedMock(real_name='test_run')
