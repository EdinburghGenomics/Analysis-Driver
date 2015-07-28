import os
import yaml
import sys

from analysis_driver.util import AppLogger


class Configuration(AppLogger):
    """
    Loads a yaml config file from the user's home, '~/.analysisdriver.yaml'
    """
    def __init__(self):
        config_file = self.__class__._find_config_file()
        self.info('Using config file at ' + config_file)
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

    def logging_config(self, log_level=None):

        dict_config = self['logging']
        dict_config['version'] = 1

        if log_level:
            assert log_level in ['DEBUG', 'INFO', 'WARN', 'ERROR', 'CRITICAL']

            for domain in ['handlers', 'loggers']:
                for k, v in dict_config[domain].items():
                    if dict_config[domain][k]['level']:
                        dict_config[domain][k]['level'] = log_level

            if dict_config['root']['level']:
                dict_config['root']['level'] = log_level

        return dict_config

x = {
    'handlers': {
        'default_hnd': {
            'formatter': 'default_fmt',
            'stream': 'ext://sys.stdout',
            'level': 'INFO'
        }
    },
    'loggers': {
        'default_logger': {
            'handlers': ['default_hnd'],
            'propagate': False,
            'level': 'INFO'
        }
    },
    'version': 1,
    'formatters': {
        'default_fmt': {
            'datefmt': '%Y-%b-%d %H:%M:%S',
            'format': '%(asctime)s [%(name)-12s][%(levelname)-5s] %(message)s'
        }
    }
}



default = Configuration()  # Singleton for access by other modules
