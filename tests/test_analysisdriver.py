__author__ = 'mwham'
from unittest import TestCase
import os.path
import json


class TestAnalysisDriver(TestCase):
    assets_path = os.path.join(os.path.dirname(__file__), 'assets')
    sample_sheet_path = os.path.join(assets_path, 'SampleSheet_analysis_driver.csv')
    barcoded_samplesheet_path = os.path.join(assets_path, 'test_runs', 'barcoded_run', 'SampleSheet_analysis_driver.csv')
    barcodeless_samplesheet_path = os.path.join(assets_path, 'test_runs', 'barcodeless_run', 'SampleSheet_analysis_driver.csv')
    barcoded_run_info_path = os.path.join(assets_path, 'test_runs', 'barcoded_run')
    barcodeless_run_info_path = os.path.join(assets_path, 'test_runs', 'barcodeless_run')
    data_output = os.path.join(assets_path, 'data_output')
    fastq_path = os.path.join(assets_path, 'fastqs')
    execs = os.path.join(assets_path, 'fake_tools')
    data_transfer = os.path.join(assets_path, 'data_transfer')

    @classmethod
    def exec_path(cls, executable):
        return os.path.join(cls.execs, executable)

    @staticmethod
    def compare_lists(observed, expected):
        if sorted(observed) != sorted(expected):
            print('')
            print('observed')
            print(observed)
            print('expected')
            print(expected)
            raise AssertionError

    def query_args_from_url(self, url):
        query_string = url.split('?')[1]
        d = {}
        for q in query_string.split('&'):
            k, v = q.split('=')
            if v.startswith('{') and v.endswith('}'):
                v = json.loads(v)
            d[k] = v

        return json.loads(json.dumps(d))
