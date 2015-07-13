import os
import yaml
from util import AnalysisDriverError


class Configuration:
    def __init__(self):
        user_config_file = os.path.join(os.getenv('HOME'), '.analysisdriver.yaml')
        local_config_file = os.path.join(os.path.dirname(__file__), 'example_analysisdriver.yaml')

        if os.path.exists(user_config_file):
            self.config_file = user_config_file
        else:
            self.config_file = local_config_file

        run_env = os.getenv('ANALYSISDRIVERENV')
        if not run_env:
            run_env = 'testing'

        self.config = yaml.load(open(self.config_file, 'r'))[run_env]

        self.work_home = self._get_param('shared', 'work_home')
        self.fastq = self._get_param('shared', 'fastq_dir')
        self.jobs = self._get_param('shared', 'jobs_dir')

        self.job_execution = self._get_param('analysisdriver', 'job_execution')
        self.bcbio = self._get_param('analysisdriver', 'bcbio')

    def _get_param(self, domain, param):
        try:
            return self.config[domain][param]
        except KeyError:
            raise AnalysisDriverError('Could not find parameter: %s/%s' % (domain, param))
