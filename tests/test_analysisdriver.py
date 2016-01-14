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

    def compare_lists(self, observed, expected):
        print()
        print('observed', '\n', observed, '\n', 'expected', '\n', expected)
        assert observed == expected
