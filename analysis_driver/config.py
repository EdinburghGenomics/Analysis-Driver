import os
import yaml
import logging
from .exceptions import AnalysisDriverError


class Configuration:
    """
    Loads a yaml config file from the user's home, '~/.analysisdriver.yaml'
    """
    def __init__(self, config_file=None):
        self._environment = None
        if not config_file:
            config_file = self.__class__._find_config_file()
        self.config_file = config_file
        try:
            self.content = yaml.safe_load(open(config_file, 'r'))[self.environment]
        except KeyError as e:
            raise AnalysisDriverError(
                'Could not load environment \'%s\' from %s' % (self.environment, self.config_file)
            ) from e
        self._validate_file_paths(self.content)

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

    def query(self, *parts, top_level=None):
        """
        Drill down into a config, e.g:
            cfg.query('logging', 'handlers', 'a_handler', 'level')
        :param str parts: Each part of the 'path' to the desired item
        :return: The relevant item if it exists in the config, else None.
        """
        if top_level is None:
            top_level = self.content
        previous_level = top_level
        item = None

        for p in parts:
            item = previous_level.get(p)
            if item:
                previous_level = item
            else:
                return None

        return item

    @property
    def environment(self):
        """
        Detects the current yaml config environment. Returns 'testing' by default. The environment is set by
        the environment variable 'ANALYSISDRIVERENV'
        """
        if not self._environment:
            self._environment = os.getenv('ANALYSISDRIVERENV', 'testing')  # default to 'testing'
        return self._environment

    def report(self):
        return yaml.safe_dump(self.content, default_flow_style=False)

    @classmethod
    def _validate_file_paths(cls, content=None):
        """
        Recursively search through the values of self.content and if the value is an absolute file path,
         assert that it exists.
        :param content: a dict, list or str (i.e. potential file path) to validate
        """
        if type(content) is dict:
            for v in content.values():
                cls._validate_file_paths(v)
        elif type(content) is list:
            for v in content:
                cls._validate_file_paths(v)
        elif type(content) is str:
            if content.startswith('/'):
                assert os.path.exists(content), 'Invalid file path: ' + content

    @staticmethod
    def _find_config_file():
        """
        Find $ANALYSISDRIVERCONFIG, ~/.analysisdriver.yaml or Analysis-Driver/etc/analysisdriver.yaml, in
        that order
        :return: Path to the config
        """
        for config in [
            os.getenv('ANALYSISDRIVERCONFIG'),
            os.path.expanduser('~/.analysisdriver.yaml'),
            os.path.join(
                os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                'etc',
                'analysisdriver.yaml'
            )
        ]:
            if config and os.path.isfile(config):
                return config
        raise AnalysisDriverError('Could not find config file in env variable, home or etc')

    @classmethod
    def _merge_dicts(cls, default, override):
        for k in set(default.keys()).union(override.keys()):
            if k in default and k in override:
                if type(default[k]) is dict and type(override[k]) is dict:
                    yield (k, dict(cls._merge_dicts(default[k], override[k])))
                else:
                    yield (k, default[k])
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


class LoggingConfiguration:
    """
    Stores logging Formatters and Handlers. self.add_handler should only be called when initially
    setting up logging.
    """
    def __init__(self):
        self.default_formatter = logging.Formatter(
            fmt=default['logging']['format'],
            datefmt=default['logging']['datefmt']
        )
        self.blank_formatter = logging.Formatter()
        self.formatter = self.default_formatter
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

    def switch_formatter(self, formatter):
        self.formatter = formatter
        for name in self.handlers:
            self.handlers[name].setFormatter(self.formatter)


# singletons for access by other modules
default = Configuration()
logging_default = LoggingConfiguration()
