import os
import yaml
from analysis_driver.util import AppLogger


class Configuration(AppLogger):
    """
    Loads a yaml config file from the user's home, '~/.analysisdriver.yaml'
    """
    def __init__(self):
        self._environment = None

        config_file = self.__class__._find_config_file()
        self.info('Using config file at ' + config_file)
        self.config_file = open(config_file, 'r')

        yaml_content = yaml.load(self.config_file)
        # Select the correct config environment
        # self.config = yaml_content[self.environment]

        self.content = yaml_content[self.environment]
        self.content['location'] = os.path.dirname(__file__)  # special case for app location
        # Merge the shared and the specific environment
        # self.content = dict(
        #     self.__class__._merge_dicts(self.config.get('shared'), self.config.get('analysisdriver'))
        # )

    def logging_config(self, debug=False):
        """
        Parse the 'logging' configuration of the yaml config to configure logging in driver.py
        :param bool debug: An option to run with logging at the 'debug' level
        :return: A dict to pass to logging.config.dictConfig
        """
        dict_config = self['logging']
        dict_config['version'] = 1

        if debug:

            for domain in ['handlers', 'loggers']:
                for k, v in dict_config[domain].items():
                    if dict_config[domain][k]['level']:
                        dict_config[domain][k]['level'] = 'DEBUG'

            if dict_config['root']['level']:
                dict_config['root']['level'] = 'DEBUG'

        return dict_config

    @property
    def environment(self):
        """
        Detects the current yaml config environment. Returns 'testing' by default. The environment is set by
        the environment variable 'ANALYSISDRIVERENV'
        """
        if not self._environment:
            self._environment = os.getenv('ANALYSISDRIVERENV', 'testing')  # Default to 'testing'
        return self._environment

    @classmethod
    def _find_config_file(cls):
        """
        Find either ~/.analysisdriver.yaml or Analysis-Driver/etc/example_analysisdriver.yaml
        :return: Path to the config
        """
        home_config = os.path.expanduser('~/.analysisdriver.yaml')
        local_config = os.path.join(os.path.dirname(__file__), 'etc', '.analysisdriver.yaml')

        if os.path.isfile(home_config):
            return home_config
        elif os.path.isfile(local_config):
            return local_config

    @classmethod
    def _merge_dicts(cls, default, override):
        """
        Recursively merges an 'override' dictionary into a 'default' dictionary.

        Adopted from http://stackoverflow.com/a/7205672/1167094

        :param default: The default dict/yaml domain
        :param override: The dict that will override values from the default
        :return: Overridden values
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
        """
        Allow access to the element of the config dict-style, e.g. config['this'] or config['this']['that']
        :param str item: A config item to retrieve
        :return: A string or deeper dict value
        """
        return self.content[item]

default = Configuration()  # singleton for access by other modules
