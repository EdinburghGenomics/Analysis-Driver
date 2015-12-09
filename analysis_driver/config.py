import os
import yaml
import logging
from .exceptions import AnalysisDriverError


class Configuration:
    def __init__(self, config_file=None):
        self.environment = os.getenv('ANALYSISDRIVERENV', 'default')
        if not config_file:
            config_file = self._find_config_file()
        self.config_file = config_file

        full_config = yaml.safe_load(open(self.config_file, 'r'))
        if not full_config.get('default'):
            raise AnalysisDriverError('Could not find \'default\' environment in ' + self.config_file)

        self.content = dict(self._merge_dicts(full_config['default'], full_config[self.environment]))

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
        Drill down into a config, e.g. cfg.query('logging', 'handlers', 'a_handler', 'level')
        :return: The relevant item if it exists in the config, else None.
        """
        if top_level is None:
            top_level = self.content
        item = None

        for p in parts:
            item = top_level.get(p)
            if item:
                top_level = item
            else:
                return None

        return item

    def report(self):
        return yaml.safe_dump(self.content, default_flow_style=False)

    @classmethod
    def validate_file_paths(cls, content=None):
        """
        Recursively search through the values of self.content and if the value is an absolute file path,
        assert that it exists.
        :param content: a dict, list or str (i.e. potential file path) to validate
        """
        invalid_file_paths = []
        if type(content) is dict:
            for v in content.values():
                invalid_file_paths.extend(cls.validate_file_paths(v))
        elif type(content) is list:
            for v in content:
                invalid_file_paths.extend(cls.validate_file_paths(v))
        elif type(content) is str:
            if content.startswith('/') and not os.path.exists(content):
                invalid_file_paths.append(content)
        return invalid_file_paths

    @staticmethod
    def _find_config_file():
        """
        Find $ANALYSISDRIVERCONFIG, ~/.analysisdriver.yaml or Analysis-Driver/etc/analysisdriver.yaml, in
        that order
        :return: Path to the config
        """
        for config in [os.getenv('ANALYSISDRIVERCONFIG'), os.path.expanduser('~/.analysisdriver.yaml')]:
            if config and os.path.isfile(config):
                return config
        raise AnalysisDriverError('Could not find config file in env variable, home or etc')

    @classmethod
    def _merge_dicts(cls, default_dict, override_dict):
        """
        Recursively merge a default dict and an overriding dict.
        """
        for k in set(override_dict.keys()).union(default_dict.keys()):
            if k in default_dict and k in override_dict:
                if type(default_dict[k]) is dict and type(override_dict[k]) is dict:
                    yield k, dict(cls._merge_dicts(default_dict[k], override_dict[k]))
                else:
                    yield k, override_dict[k]
            elif k in default_dict:
                yield k, default_dict[k]
            else:
                yield k, override_dict[k]

    def merge(self, override_dict):
        """
        Merge the provided dict with the config content potententially overiding existing parameters
        """
        self.content = dict(self._merge_dicts(self.content, override_dict))

    def __getitem__(self, item):
        """
        Allow dict-style access, e.g. config['this'] or config['this']['that']
        """
        return self.content[item]

    def __contains__(self, item):
        """
        Allow search in the first layer of the config with "in" operator
        """
        return self.content.__contains__(item)


class LoggingConfiguration:
    """
    Stores logging Formatters and Handlers.
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
        :param str name: A name or id to assign the Handler
        :param logging.FileHandler handler:
        """
        handler.setFormatter(self.formatter)
        handler.setLevel(self.log_level)
        self.handlers[name] = handler

    def switch_formatter(self, formatter):
        """
        Set all handlers to formatter
        """
        self.formatter = formatter
        for name in self.handlers:
            self.handlers[name].setFormatter(self.formatter)


# singletons for access by other modules
default = Configuration()
logging_default = LoggingConfiguration()
