__author__ = 'mwham'
from tests.test_analysisdriver import TestAnalysisDriver
from analysis_driver import executor


class TestExecutor(TestAnalysisDriver):
    def setUp(self):
        self.executor = executor.Executor(['ls', self.assets_path])

    def test_cmd(self):
        assert self.executor.cmd == ['ls', self.assets_path]

    def test_process(self):
        pass