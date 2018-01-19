from os.path import join, dirname
from unittest import TestCase
from unittest.mock import Mock
from egcg_core.app_logging import logging_default as log_cfg
from analysis_driver.config import default as cfg, etc_config
from analysis_driver import tool_versioning


class NamedMock(Mock):
    @property
    def name(self):
        return self.real_name


class TestAnalysisDriver(TestCase):
    assets_path = join(dirname(__file__), 'assets')
    fastq_path = join(assets_path, 'fastqs')


cfg.load_config_file(etc_config('example_analysisdriver.yaml'))
helper = TestAnalysisDriver()
log_cfg.cfg = cfg['logging']
log_cfg.add_stdout_handler()

tool_versioning.toolset = Mock(__getitem__=lambda self, key: 'path/to/' + key)
