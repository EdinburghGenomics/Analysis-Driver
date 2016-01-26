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
        print()
        print('observed', '\n', observed, '\n', 'expected', '\n', expected)
        assert sorted(observed) == sorted(expected)

    def setup_db(self, db, endpoints):
        db.clear()
        for e in endpoints:
            db[e] = self._fake_rest_data(e)

    def _fake_rest_data(self, endpoint):
        return {
            'data': [],
            '_meta': {'max_results': 25, 'total': 0, 'page': 1},
            '_links': {
                'self': {'href': endpoint, 'title': endpoint},
                'parent': {'href': '/', 'title': 'home'}
            }
        }
