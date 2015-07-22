import os
import yaml
import sys

from analysis_driver.util import AppLogger, AnalysisDriverError


class Configuration(AppLogger):
    """
    Loads a yaml config file from the user's home, '~/.analysisdriver.yaml'
    """
    def __init__(self):
        config_file = self.__class__._find_config_file()
        self.log('Using config file at ' + config_file)
        self.config_file = open(config_file, 'r')

        yaml_content = yaml.load(self.config_file)
        self._environment = None
        # Select the correct config environment
        self.config = yaml_content[self.environment]

        # Merge the shared and the specific environment
        self.content = dict(
            self.__class__._merge_dicts(self.config.get('shared'), self.config.get('analysisdriver'))
        )

    @property
    def environment(self):
        """
        Detects the current environment. Returns 'testing' by default.
        """
        if not self._environment:
            self._environment = os.getenv('ANALYSISDRIVERENV', 'testing')  # Default to 'testing'
        return self._environment

    @classmethod
    def _find_config_file(cls):
        home_config = os.path.expanduser('~/.analysisdriver.yaml')
        local_config = os.path.join(os.path.dirname(__file__), 'example_analysisdriver.yaml')

        if os.path.isfile(home_config):
            return home_config
        elif os.path.isfile(local_config):
            return local_config

    @classmethod
    def _merge_dicts(cls, default, override):
        """
        Recursively merges an 'override' dictionary into a 'default' dictionary.

        Adopted from http://stackoverflow.com/a/7205672/1167094
        """
        for k in set(default.keys()).union(override.keys()):
            if k in default and k in override:
                if isinstance(default[k], dict) and isinstance(override[k], dict):
                    yield (k, dict(cls._merge_dicts(default[k], override[k])))
                else:
                    yield (k, override[k])
            elif k in default:
                yield (k, default[k])
            else:
                yield (k, override[k])

    def __getitem__(self, item):
        # Allow access to the element of the config with dictionary style
        return self.content[item]

    def logging(self, log_level='INFO'):
        logging_config = self['logging']
        dict_config = {
            'version': 1,
            'formatters': {},
            'handlers': {},
            'loggers': {
                '': {
                    'handlers': [],
                    'propagate': False
                }
            }
        }
        for name, formatter in logging_config['formatters'].items():
            dict_config['formatters'][name] = formatter

        for name, handler in logging_config['handlers'].items():
            dict_config['handlers'][name] = {
                'level': handler['level'],
                'formatter': handler['formatter'],
                'stream': handler['stream'],
                'class': 'logging.StreamHandler'
            }
            if 'sys.stdout' in handler['stream']:  # Special case if the user wants to log to sys.stdout
                dict_config['handlers'][name]['stream'] = sys.stdout

            dict_config['loggers']['']['handlers'].append(name)

        dict_config['loggers']['']['level'] = logging_config['default_level']
        return dict_config


from unittest import TestCase


class TestConfig(TestCase):
    def setUp(self):
        self.config = Configuration()

    def test_getitem(self):
        self.assertEqual(self.config['raw_dir'], 'raw')
