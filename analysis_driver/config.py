import os
import yaml
from .app_logging import AppLogger
from .exceptions import AnalysisDriverError


class Configuration(AppLogger):
    """
    Loads a yaml config file from the user's home, '~/.analysisdriver.yaml'
    """
    def __init__(self):
        self._environment = None

        config_file = self.__class__._find_config_file()
        self.info('Using config file at ' + config_file)
        self.config_file = open(config_file, 'r')

        try:
            self.content = yaml.load(self.config_file)[self.environment]
        except KeyError as e:
            raise AnalysisDriverError('Could not load environment \'%s\'' % self.environment) from e
        self.content['location'] = os.path.dirname(__file__)  # special case for app location

    def logging_config(self, debug=False, d_handler=None, no_stdout=False):
        """
        Parse the 'logging' configuration of the yaml config to configure logging in driver.py
        :param bool debug: An option to run with logging at the 'debug' level
        :return: A dict to pass to logging.config.dictConfig
        """
        dict_config = self['logging']
        dict_config['version'] = 1

        if d_handler:
            dict_config['handlers']['d_handler'] = d_handler
            dict_config['root']['handlers'].append('d_handler')

        if no_stdout:
            for name, handler in dict_config['handlers'].items():
                if 'stream' in handler:
                    if handler['stream'] == 'ext://sys.stdout' or handler['stream'] == 'ext://sys.stderr':
                        handler['stream'] = open(os.devnull, 'w')

        if debug:
            for domain in ['handlers', 'loggers']:
                if domain in dict_config:
                    for k, v in dict_config[domain].items():
                        if dict_config[domain][k]['level']:
                            dict_config[domain][k]['level'] = 'DEBUG'
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

    def __getitem__(self, item):
        """
        Allow access to the element of the config dict-style, e.g. config['this'] or config['this']['that']
        :param str item: A config item to retrieve
        :return: A string or deeper dict value
        """
        return self.content[item]

default = Configuration()  # singleton for access by other modules
