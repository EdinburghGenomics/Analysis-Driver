from unittest import TestCase
from unittest.mock import Mock

from egcg_core.app_logging import logging_default as log_cfg
from analysis_driver.config import default as cfg, _etc_config
from os.path import join, dirname
import json

class NamedMock(Mock):

    @property
    def name(self):
        return self.real_name


class TestAnalysisDriver(TestCase):
    assets_path = join(dirname(__file__), 'assets')

    data_output = join(assets_path, 'data_output')
    fastq_path = join(assets_path, 'fastqs')
    execs = join(assets_path, 'fake_tools')
    data_transfer = join(assets_path, 'data_transfer')

    def __init__(self, *args, **kwargs):
        super(TestAnalysisDriver, self).__init__(*args, **kwargs)
        cfg.load_config_file(_etc_config('example_analysisdriver.yaml'))

    @classmethod
    def exec_path(cls, executable):
        return join(cls.execs, executable)

    @staticmethod
    def compare_lists(observed, expected):
        if sorted(observed) != sorted(expected):
            print('')
            print('observed')
            print(observed)
            print('expected')
            print(expected)
            raise AssertionError

    @staticmethod
    def query_args_from_url(url):
        query_string = url.split('?')[1]
        d = {}
        for q in query_string.split('&'):
            k, v = q.split('=')
            if v.startswith('{') and v.endswith('}'):
                v = json.loads(v)
            d[k] = v

        return json.loads(json.dumps(d))


helper = TestAnalysisDriver()
log_cfg.cfg = cfg['logging']
log_cfg.add_stdout_handler()
