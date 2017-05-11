import luigi
import json
from os.path import join
from egcg_core.app_logging import AppLogger
from analysis_driver.exceptions import PipelineError
from analysis_driver.config import default as cfg


class EGCGParameter(luigi.Parameter):
    """Parameter that does not call `warnings.warn` in self.serialize."""
    def serialize(self, x):
        if type(x) in (list, tuple, dict, set):
            x = '<%s of len %s>' % (x.__class__.__name__, len(x))
        return str(x)

    # def __getattribute__(self, item):
    #     return getattr(self, item, None)


class EGCGListParameter(luigi.ListParameter):
    def serialize(self, x):
        return json.dumps(x, cls=EGCGEncoder)


class EGCGEncoder(json.JSONEncoder):
    def default(self, o):
        return str(o)


class BasicStage(luigi.Task, AppLogger):
    __stagename__ = None
    exit_status = None
 
    previous_stages = EGCGListParameter(default=[])
    dataset = EGCGParameter()

    def output(self):
        return [luigi.LocalTarget(join(self.job_dir, '.' + self.stage_name + '.stage'))]

    @property
    def stage_name(self):
        return self.__stagename__ or self.__class__.__name__.lower()

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


class Stage(BasicStage):
    def output(self):
        return [RestAPITarget(self)]

    def run(self):
        self.dataset.start_stage(self.stage_name)
        self.exit_status = self._run()
        self.dataset.end_stage(self.stage_name, self.exit_status)
        if self.exit_status:
            raise PipelineError('Exit status was %s. Stopping' % self.exit_status)

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


# example Luigi workflow
from os.path import dirname
from tests.test_analysisdriver import NamedMock


class Stage1(Stage):
    def _run(self):
        print('Doing some things')
        return 0


class Stage2(Stage1):
    pass


class Stage25(Stage1):
    __stagename__ = 'stage2'


class Stage3(Stage1):
    pass


def example_pipeline(d):
    stage1 = Stage1(dataset=d)
    stage2 = Stage2(dataset=d, previous_stages=[stage1])
    stage2_5 = Stage25(dataset=d, previous_stages=[stage1])
    stage3 = Stage3(dataset=d, previous_stages=[stage1, stage2, stage2_5])

    stage3.exit_status = 9
    luigi.build([stage3], local_scheduler=True)
    return stage3.exit_status


if __name__ == '__main__':
    cfg.load_config_file(join(dirname(dirname(__file__)), 'etc', 'example_analysisdriver.yaml'))
    dataset = NamedMock(real_name='test_dataset')
    print(example_pipeline(dataset))
