__author__ = 'mwham'
from unittest import TestCase
import os.path


class TestAnalysisDriver(TestCase):
    assets_path = os.path.join(os.path.dirname(__file__), 'assets')
    sample_sheet_path = os.path.join(assets_path, 'SampleSheet_analysis_driver.csv')
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
