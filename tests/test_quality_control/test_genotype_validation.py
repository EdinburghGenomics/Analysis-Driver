from tests.test_analysisdriver import TestAnalysisDriver

__author__ = 'tcezard'


class TestGenotypeValidation(TestAnalysisDriver):
    def setUp(self):
        self.executor = executor.Executor(['ls', self.assets_path])

    def test_cmd(self):
        assert self.executor.cmd == ['ls', self.assets_path]

    def test_process(self):
        pass