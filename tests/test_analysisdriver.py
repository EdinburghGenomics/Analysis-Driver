from unittest import TestCase
from egcg_core.app_logging import logging_default as log_cfg
from analysis_driver.config import default as cfg, load_config, _etc_config
from os.path import join, dirname
import json


class TestAnalysisDriver(TestCase):
    assets_path = join(dirname(__file__), 'assets')
    sample_sheet_path = join(assets_path, 'SampleSheet_analysis_driver.csv')
    barcoded_samplesheet_path = join(assets_path, 'test_runs', 'barcoded_run', 'SampleSheet_analysis_driver.csv')
    barcodeless_samplesheet_path = join(assets_path, 'test_runs', 'barcodeless_run', 'SampleSheet_analysis_driver.csv')
    barcoded_run_info_path = join(assets_path, 'test_runs', 'barcoded_run')
    barcodeless_run_info_path = join(assets_path, 'test_runs', 'barcodeless_run')
    data_output = join(assets_path, 'data_output')
    fastq_path = join(assets_path, 'fastqs')
    execs = join(assets_path, 'fake_tools')
    data_transfer = join(assets_path, 'data_transfer')

    def __init__(self, *args, **kwargs):
        super(TestAnalysisDriver, self).__init__(*args, **kwargs)
        load_config(_etc_config('example_analysisdriver.yaml'))

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
