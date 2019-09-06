import os
import random
import shutil
import string

import luigi
import json
from os.path import join
from egcg_core.app_logging import AppLogger
from analysis_driver.exceptions import PipelineError
from analysis_driver.config import default as cfg


class Parameter(luigi.Parameter):
    """Parameter that does not call `warnings.warn` in self.serialize."""
    def serialize(self, x):
        if type(x) in (list, tuple, dict, set):
            x = '<%s of len %s>' % (x.__class__.__name__, len(x))
        return str(x)

    # def __getattribute__(self, item):
    #     return getattr(self, item, None)


class ListParameter(luigi.ListParameter):
    def serialize(self, x):
        return json.dumps(x, cls=Encoder)


class Encoder(json.JSONEncoder):
    def default(self, o):
        return str(o)


class BasicStage(luigi.Task, AppLogger):
    exit_status = None
    retry_count = 1  # turn off automatic retrying upon task failure

    stage_name = Parameter(default=None)
    previous_stages = ListParameter(default=[])
    dataset = Parameter()
    pipeline = Parameter(default=None)

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.stage_name = self.stage_name or self.__class__.__name__.lower()

    def output(self):
        return [luigi.LocalTarget(join(self.job_dir, '.' + self.stage_name + '.stage'))]

    def requires(self):
        for s in self.previous_stages:
            yield s

    @property
    def job_dir(self):
        return join(cfg['jobs_dir'], self.dataset.name)

    def __str__(self):
        return '%s(previous_stages=%s, dataset=%s)' % (
            self.__class__.__name__, [s.__class__.__name__ for s in self.previous_stages], self.dataset
        )

    def on_failure(self, exception):
        self.dataset.register_exception(self, exception)
        if self.exit_status is None:
            self.dataset.end_stage(self.stage_name, 9)
        return super().on_failure(exception)


class Stage(BasicStage):

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.dir_to_delete = []

    def output(self):
        return [RestAPITarget(self)]

    def run(self):
        self.dataset.start_stage(self.stage_name)
        self.exit_status = self._run()
        self.dataset.end_stage(self.stage_name, self.exit_status)
        if self.exit_status:
            raise PipelineError('Exit status for %s was %s. Stopping' % (self.stage_name, self.exit_status))
        # Remove temp directories
        for d in self.dir_to_delete:
            shutil.rmtree(d)
        self.info('Finished stage %s', self.stage_name)

    @property
    def input_dir(self):
        return join(cfg[self.dataset.type]['input_dir'], self.dataset.name)

    def _run(self):
        return 0

    @property
    def exec_dir(self):
        """Directory where the slurm jobs and slurm logs are located."""
        d = os.path.join(self.job_dir, 'slurm_and_logs')
        os.makedirs(d, exist_ok=True)
        return d

    @staticmethod
    def _id_generator(size=6, chars=string.ascii_lowercase + string.digits):
        return ''.join(random.choice(chars) for _ in range(size))

    def create_tmp_dir(self):
        """Create a temporary directory that will be used during this stage and deleted at the end."""
        tmp_dir = os.path.join(self.job_dir, self._id_generator())
        os.makedirs(tmp_dir)
        self.dir_to_delete.append(tmp_dir)
        return tmp_dir


class RestAPITarget(luigi.Target):
    def __init__(self, stage):
        self.stage = stage

    @property
    def rest_api_stage(self):
        return self.stage.dataset.get_stage(self.stage.stage_name)

    def exists(self):
        s = self.rest_api_stage
        return s and bool(s.get('date_finished')) and s.get('exit_status') == 0
