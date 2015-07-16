import logging
import os
import yaml


home = os.getenv('HOME')
user_config_file = os.path.join(home, '.analysisdriver.yaml')
local_config_file = os.path.join(os.path.dirname(__file__), 'example_analysisdriver.yaml')

if os.path.exists(user_config_file):
    config_file = user_config_file
else:
    config_file = local_config_file

run_env = os.getenv('ANALYSISDRIVERENV', 'testing')

config = yaml.load(open(config_file, 'r'))[run_env]

work_home = config['shared']['work_home']
fastq = config['shared']['fastq_dir']
jobs = config['shared']['jobs_dir']

job_execution = config['analysisdriver']['job_execution']


class Config():
    """
    Loads a configuration file and returns a Config object
    """
    def __init__(self):
        self.config_file = open(self.__class__._find_config_file(), 'r')
        tmp_config = yaml.load(self.config_file)
        self._environment = None
        #Select the config specific to the environment
        self.config = tmp_config[self.environment]

        #Merge the shared and the specific environment
        self.content = dict(
                self.__class__._merge_dicts(self.config.get('shared'), self.config.get('analysisdriver'))
            )

    @property
    def environment(self):
        """
        Detects the current environment. Returns 'testing' by default.
        """
        if not self._environment:
            self._environment = os.getenv('ANALYSISDRIVERENV', 'testing')
        return self._environment

    @classmethod
    def _find_config_file(cls):
        possible_paths = [os.path.join(home, '.analysisdriver.yaml'),
                          os.path.join(os.path.dirname(__file__), 'example_analysisdriver.yaml')]
        for p in possible_paths:
            if p and os.path.isfile(p):
                logging.info('Found config file at {}'.format(p))
                return p

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
        #Allow access to the element of the config with dictionary style
        return self.content[item]


    def __getattr__(self, item):
        #Allow access to the first layer with attribute style
        return self.content[item]



from unittest import TestCase

class test_config(TestCase):
    def setUp(self):
        self.config=Config()

    def test_getitem(self):
        self.assertEqual(self.config['raw_dir'],'raw')

    def test_getattr(self):
        self.assertEqual(self.config.raw_dir,'raw')