__author__ = 'mwham'
import os
from tests.test_analysisdriver import TestAnalysisDriver
from analysis_driver.config import Configuration, EnvConfiguration, LoggingConfiguration, logging_default
from analysis_driver import app_logging
from analysis_driver.exceptions import AnalysisDriverError
import pytest
import logging
import sys


class TestConfiguration(TestAnalysisDriver):
    def setUp(self):
        self.cfg = EnvConfiguration(
            (
                os.getenv('ANALYSISDRIVERCONFIG'),
                os.path.join(self.assets_path, '..', '..', 'etc', 'example_analysisdriver.yaml')
            )
        )
        self.sample_sheet_cfg = Configuration(
            [os.path.join(self.assets_path, '..', '..', 'etc', 'sample_sheet_cfg.yaml')]
        )

    def test_get(self):
        get = self.cfg.get
        assert get('tt_agent_delay') == 20

        assert get('nonexistent_thing') is None
        assert get('nonexistent_thing', 'a_default') == 'a_default'
        assert self.cfg.get('logging').get('handlers').get('stdout').get('level') == 'INFO'

    def test_find_config_file(self):
        existing_cfg_file = os.path.join(self.assets_path, '..', '..', 'etc', 'example_analysisdriver.yaml')
        non_existing_cfg_file = os.path.join(self.assets_path, 'a_file_that_does_not_exist.txt')
        assert self.cfg._find_config_file((non_existing_cfg_file, existing_cfg_file)) == existing_cfg_file

        with pytest.raises(AnalysisDriverError) as e:
            self.cfg._find_config_file((non_existing_cfg_file,))
            assert 'Could not find config file in self.cfg_search_path' in str(e)

    def test_query(self):
        assert self.cfg.query('logging', 'handlers', 'stdout', 'level') == 'INFO'
        assert self.cfg.query('nonexistent_thing') is None
        assert self.cfg.query('logging', 'handlers', 'nonexistent_handler') is None
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
            self.sample_sheet_cfg.query(
                'sample_project',
                top_level=self.sample_sheet_cfg.query('column_names')
            ),
            self.sample_sheet_cfg.query(
                'column_names',
                'sample_project'
            ),
            self.sample_sheet_cfg.query(
                'sample_project',
                top_level=self.sample_sheet_cfg.query('column_names')
            )
        ]:
            assert q == ['Sample_Project', 'SampleProject', 'Project_Name']

    def test_merge_dicts(self):
        default_dict = {
            'this': {
                'another': [2, 3, 4],
                'more': {
                    'thing': 'thang'
                }
            },
            'that': 'a_thing',
            'other': {
                'another': [2, '3', 4],
                'more': {
                    'thing': 'thang'
                }
            }
        }
        override_dict = {
            'that': 'another_thing',
            'another': 4,
            'more': 5,
            'other': {
                'another': [8, 9, 10],
                'more': {'thung': 'theng'}
            }
        }
        merged_dict = self.cfg._merge_dicts(default_dict, override_dict)

        assert dict(merged_dict) == {
            'this': {
                'another': [2, 3, 4],
                'more': {
                    'thing': 'thang'
                }
            },
            'that': 'another_thing',
            'other': {
                'another': [8, 9, 10],
                'more': {
                    'thing': 'thang',
                    'thung': 'theng'
                }
            },
            'another': 4,
            'more': 5
        }

        assert dict(self.cfg._merge_dicts(default_dict, default_dict)) == default_dict


class TestLoggingConfiguration(TestAnalysisDriver):
    def setUp(self):
        self.log_cfg = LoggingConfiguration()

    def test_formatters(self):
        default = self.log_cfg.default_formatter
        assert default.datefmt == '%Y-%b-%d %H:%M:%S'
        assert default._fmt == '[%(asctime)s][%(name)s][%(levelname)s] %(message)s'

        blank = self.log_cfg.blank_formatter
        assert blank.datefmt is None
        assert blank._fmt == '%(message)s'

        assert self.log_cfg.formatter is default

        assert self.log_cfg.handlers == {}
        assert self.log_cfg.default_level == logging.INFO

    def test_add_handler(self):
        h = logging.StreamHandler(stream=sys.stdout)
        self.log_cfg.add_handler('test', h)

        assert self.log_cfg.handlers['test'] is h
        assert h.formatter is self.log_cfg.formatter
        assert h.level == logging.INFO

    def test_switch_formatter(self):
        h = logging.StreamHandler(stream=sys.stdout)
        self.log_cfg.add_handler('test', h)

        self.log_cfg.switch_formatter(self.log_cfg.blank_formatter)
        assert h.formatter is self.log_cfg.blank_formatter
        self.log_cfg.switch_formatter(self.log_cfg.default_formatter)
        assert h.formatter is self.log_cfg.default_formatter


class TestAppLogging(app_logging.AppLogger, TestAnalysisDriver):
    def setUp(self):
        logging_default.add_handler('test', logging.StreamHandler(stream=sys.stdout))
        logging_default.add_handler('test2', logging.StreamHandler(stream=sys.stderr))

    def tearDown(self):
        logging_default.handlers = {}

    def test_log_msgs(self):
        self.debug('Debug')
        self.info('Info')
        self.warn('Warning')
        self.error('Error')
        self.critical('Critical')

        with pytest.raises(KeyError) as e:
            self.critical('Oh noes!', error_class=KeyError)
            assert str(e) == 'Oh noes!'

    def test_get_logger(self):
        logger = app_logging.get_logger('test')
        assert logger.level == logging_default.default_level
        assert list(logging_default.handlers.values()) == logger.handlers
