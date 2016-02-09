import os
import yaml
import logging
from .exceptions import AnalysisDriverError


class Configuration:
    def __init__(self, cfg_search_path):
        self.cfg_search_path = cfg_search_path
        self.config_file = self._find_config_file()
        self.content = yaml.safe_load(open(self.config_file, 'r'))

    def _find_config_file(self):
        for p in self.cfg_search_path:
            if p and os.path.isfile(p):
                return p
        raise AnalysisDriverError('Could not find config file in self.cfg_search_path')

    def get(self, item, ret_default=None):
        """
        Dict-style item retrieval with default
        :param item: The key to search for
        :param ret_default: What to return if the key is not present
        """
        try:
            return self[item]
        except KeyError:
            return ret_default

    def query(self, *parts, top_level=None, ret_default=None):
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
                return ret_default

        return item

    def report(self):
        return yaml.safe_dump(self.content, default_flow_style=False)

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


class EnvConfiguration(Configuration):
    def __init__(self, cfg_search_path):
        super().__init__(cfg_search_path)
        env = os.getenv('ANALYSISDRIVERENV', 'default')
        if not self.content.get('default'):
            raise AnalysisDriverError('Could not find \'default\' environment in ' + self.config_file)
        self.content = dict(self._merge_dicts(self.content['default'], self.content[env]))

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


class LoggingConfiguration:
    """
    Stores logging Formatters and Handlers.
    """
    def __init__(self):
        self.default_formatter = logging.Formatter(
            fmt=default.query(
                'logging',
                'format',
                ret_default='[%(asctime)s][%(name)s][%(levelname)s] %(message)s'
            ),
            datefmt=default.query('logging', 'datefmt', ret_default='%Y-%b-%d %H:%M:%S')
        )
        self.blank_formatter = logging.Formatter()
        self.formatter = self.default_formatter
        self.handlers = {}
        self.default_level = logging.INFO

    def add_handler(self, name, handler, level=None):
        """
        :param str name: A name or id to assign the Handler
        :param logging.FileHandler handler:
        """
        if level is None:
            level = self.default_level
        handler.setFormatter(self.formatter)
        handler.setLevel(level)
        self.handlers[name] = handler

    def switch_formatter(self, formatter):
        """
        Set all handlers to formatter
        """
        self.formatter = formatter
        for name in self.handlers:
            self.handlers[name].setFormatter(self.formatter)


def _dir_path():
    """Find the absolute path of 2 dirs above this file (should be Analysis-Driver)"""
    return os.path.dirname(os.path.abspath(os.path.dirname(__file__)))


def _etc_config(config_file):
    return os.path.join(_dir_path(), 'etc', config_file)


# singletons for access by other modules
default = EnvConfiguration(
    [
        os.getenv('ANALYSISDRIVERCONFIG'),
        os.path.expanduser('~/.analysisdriver.yaml'),
        os.path.join(os.path.dirname(os.path.dirname(__file__)), 'etc', 'example_analysisdriver.yaml')

    ]
)
output_files_config = Configuration([_etc_config('output_files.yaml')])
sample_sheet_config = Configuration([_etc_config('sample_sheet_cfg.yaml')])
logging_default = LoggingConfiguration()
