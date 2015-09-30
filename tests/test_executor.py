__author__ = 'mwham'
from tests.test_analysisdriver import TestAnalysisDriver
from analysis_driver import executor
import os.path
import sys
import logging

logging.basicConfig(stream=sys.stdout, level=logging.DEBUG)


class TestExecutor(TestAnalysisDriver):
    def setUp(self):
        self.executor = executor.Executor(['ls', self.assets_path])

    def test_cmd(self):
        assert self.executor.cmd == ['ls', self.assets_path]

    def test_stream(self):
        e = executor.StreamExecutor(
            [os.path.join(os.path.dirname(__file__), 'assets', 'countdown.sh')]
        )
        e.start()
        assert e.join() == 0

        f = executor.StreamExecutor(
            [os.path.join(os.path.dirname(__file__), 'assets', 'countdown.sh'), 'dodgy'],
        )
        f.start()
        assert f.join() == 13
