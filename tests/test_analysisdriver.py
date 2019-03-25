from os.path import join, dirname, abspath
from unittest import TestCase
from unittest.mock import Mock
from analysis_driver.config import default as cfg, etc_config
from analysis_driver import tool_versioning


class NamedMock(Mock):
    @property
    def name(self):
        return self.real_name


class TestAnalysisDriver(TestCase):
    assets_path = join(abspath(dirname(__file__)), 'assets')
    fastq_path = join(assets_path, 'fastqs')

    @staticmethod
    def _touch(input_file):
        open(input_file, 'w').close()

    @classmethod
    def setUpClass(cls):
        cfg.load_config_file(etc_config('example_analysisdriver.yaml'))
        tool_versioning.toolset.versioning_cfg.load_config_file(join(cls.assets_path, 'tool_versioning.yaml'))
        tool_versioning.toolset.select_type('fake_tools')
        tool_versioning.toolset.select_version(0)
