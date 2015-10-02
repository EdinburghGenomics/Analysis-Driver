__author__ = 'mwham'
from tests.test_analysisdriver import TestAnalysisDriver
from analysis_driver import executor
import os.path
import sys
import logging

logging.basicConfig(stream=sys.stdout, level=logging.DEBUG)


class TestExecutor(TestAnalysisDriver):
    def test_cmd(self):
        e = executor.Executor(['ls', os.path.join(self.assets_path, '..')])
        out, err = e.run()
        for f in ['test_executor.py', 'test_util.py', 'test_notification', '__init__.py']:
            assert f in str(out)

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
