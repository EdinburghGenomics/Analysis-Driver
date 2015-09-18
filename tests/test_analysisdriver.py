__author__ = 'mwham'
from unittest import TestCase
import os.path
from analysis_driver.config import Configuration, LoggingConfiguration


class TestAnalysisDriver(TestCase):
    assets_path = os.path.join(os.path.dirname(__file__), 'assets')
    data_output = os.path.join(assets_path, 'data_output')
    fastq_path = os.path.join(assets_path, 'fastqs')


class TestConfiguration(TestAnalysisDriver):
    def setUp(self):
        self.cfg = Configuration()

    def test_get(self):
        get = self.cfg.get
        assert get('tt_agent_delay') == 20

        assert get('nonexistent_thing') is None
        assert get('nonexistent_thing', 'a_default') == 'a_default'

        # test Configuration.get's compatibility with dict.get
        assert self.cfg.get('logging').get('handlers').get('stdout').get('level') == 'INFO'

    def test_query(self):
        assert self.cfg.query('logging', 'handlers', 'stdout', 'level') == 'INFO'

        assert self.cfg.query('nonexistent_thing') is None
        assert self.cfg.query('logging') == {
            'handlers': {
                'stdout': {
                    'stream': 'ext://sys.stdout',
                    'level': 'INFO'
                }
            },
            'datefmt': '%Y-%b-%d %H:%M:%S',
            'format': '[%(asctime)s][%(name)s][%(levelname)s] %(message)s'
        }

        for q in [
            self.cfg.query(
                'sample_project',
                top_level=self.cfg.query('sample_sheet', 'column_names')
            ),
            self.cfg.query(
                'column_names',
                'sample_project',
                top_level=self.cfg.query('sample_sheet')
            ),
            self.cfg.query(
                'sample_project',
                top_level=self.cfg.query('column_names', top_level=self.cfg.query('sample_sheet'))
            )
        ]:
            assert q == ['Sample_Project', 'SampleProject', 'Project_Name']

