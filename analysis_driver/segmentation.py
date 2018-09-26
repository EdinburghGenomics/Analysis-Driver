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
    def output(self):
        return [RestAPITarget(self)]

    def run(self):
        self.dataset.start_stage(self.stage_name)
        self.exit_status = self._run()
        self.dataset.end_stage(self.stage_name, self.exit_status)
        if self.exit_status:
            raise PipelineError('Exit status for %s was %s. Stopping' % (self.stage_name, self.exit_status))

        self.info('Finished stage %s' % self.stage_name)

    @property
    def input_dir(self):
        return join(cfg['input_dir'], self.dataset.name)

    def _run(self):
        return 0


class RestAPITarget(luigi.Target):
    def __init__(self, stage):
        self.stage = stage

    @property
    def rest_api_stage(self):
        return self.stage.dataset.get_stage(self.stage.stage_name)

    def exists(self):
        s = self.rest_api_stage
        return s and bool(s.get('date_finished')) and s.get('exit_status') == 0
