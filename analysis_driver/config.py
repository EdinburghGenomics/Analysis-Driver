import os
import yaml
import logging
from .exceptions import AnalysisDriverError


class Configuration:
    """
    Loads a yaml config file from the user's home, '~/.analysisdriver.yaml'
    """
    def __init__(self):
        self._environment = None
        config_file = self.__class__._find_config_file()
        self.config_file = open(config_file, 'r')

        try:
            self.content = yaml.load(self.config_file)[self.environment]
        except KeyError as e:
            raise AnalysisDriverError('Could not load environment \'%s\'' % self.environment) from e

    def get(self, item, return_default=None):
        """
        Dict-style item retrieval with default
        :param item: The key to search for
        :param return_default: What to return if the key is not present
        """
        try:
            return self[item]
        except KeyError:
            return return_default

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
        local_config = os.path.join(os.path.dirname(os.path.dirname(__file__)), 'etc', '.analysisdriver.yaml')
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


class LoggingConfiguration:
    """
    Stores logging Formatters and Handlers. self.add_handler should only be called when initially
    setting up logging.
    """
    def __init__(self):
        self.formatter = logging.Formatter(
            fmt=default['logging']['format'],
            datefmt=default['logging']['datefmt']
        )
        self.handlers = {}
        self.log_level = logging.INFO

    def add_handler(self, name, handler):
        """
        :param str name:
        :param logging.FileHandler handler:
        """
        handler.setFormatter(self.formatter)
        handler.setLevel(self.log_level)
        self.handlers[name] = handler

default = Configuration()  # singleton for access by other modules
logging_default = LoggingConfiguration()
