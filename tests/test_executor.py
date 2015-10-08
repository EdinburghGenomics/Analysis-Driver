__author__ = 'mwham'
from tests.test_analysisdriver import TestAnalysisDriver
from analysis_driver import executor
from analysis_driver.exceptions import AnalysisDriverError
import pytest
import os.path
import sys
import logging

logging.basicConfig(stream=sys.stdout, level=logging.DEBUG)


class TestExecutor(TestAnalysisDriver):
    def test_cmd(self):
        e = executor.SimpleExecutor('ls ' + os.path.join(self.assets_path, '..'))
        exit_status = e.join()
        assert exit_status == 0

    def test_stream(self):
        e = executor.StreamExecutor(
            os.path.join(os.path.dirname(__file__), 'assets', 'countdown.sh')
        )
        e.start()
        assert e.join() == 0

        f = executor.StreamExecutor(
            os.path.join(os.path.dirname(__file__), 'assets', 'countdown.sh') + ' dodgy',
        )
        f.start()
        assert f.join() == 13

    def test_dodgy_stream(self):
        with pytest.raises(AnalysisDriverError) as err:
            e = executor.SimpleExecutor('dodgy_cmd')
            e.join()
            assert 'Command failed: \'dodgy_cmd\'' in str(err)

        with pytest.raises(AnalysisDriverError) as err2:
            f = executor.StreamExecutor('dodgy_cmd')
            f.start()
            f.join()
            assert 'self.proc command failed: \'dodgy_cmd\'' in str(err2)
