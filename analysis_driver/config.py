import os
import yaml
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
        self.content['location'] = os.path.dirname(__file__)  # special case for app location

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

default = Configuration()  # singleton for access by other modules
