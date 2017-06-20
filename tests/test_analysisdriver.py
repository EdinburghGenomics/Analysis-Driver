import json
from os.path import join, dirname
from unittest import TestCase
from unittest.mock import Mock
from egcg_core.app_logging import logging_default as log_cfg
from analysis_driver.config import default as cfg, _etc_config


class NamedMock(Mock):
    @property
    def name(self):
        return self.real_name


class TestAnalysisDriver(TestCase):
    assets_path = join(dirname(__file__), 'assets')
    fastq_path = join(assets_path, 'fastqs')

    def __init__(self, *args, **kwargs):
        super(TestAnalysisDriver, self).__init__(*args, **kwargs)
        cfg.load_config_file(_etc_config('example_analysisdriver.yaml'))

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
