import luigi
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


class BasicStage(luigi.Task, AppLogger):
    __stagename__ = None
    exit_status = None

    previous_stages = luigi.ListParameter(default=[])
    dataset = EGCGParameter()

    def output(self):
        return [luigi.LocalTarget(join(self.job_dir, '.' + self.stage_name + '.stage'))]

    @property
    def stage_name(self):
        return self.__stagename__ or self.__class__.__name__.lower()

    def requires(self):
        """
        Generates prior Stages from self.previous_stages, which should be a list of Stage configs:
        [
            Stage1,
            (Stage2, {'a_parameter': EGCGParameter()})
        ]
        Stage1 will have no extra params, Stage2 will be created with the params in the second tuple element.
        All Stages will be created with dataset=self.dataset
        """
        for s in self.previous_stages:
            if isinstance(s, type):
                yield s(dataset=self.dataset)
            elif type(s) in (tuple, list):
                cls, config = s
                yield cls(dataset=self.dataset, **config)

    @property
    def job_dir(self):
        return join(cfg['jobs_dir'], self.dataset.name)


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
from analysis_driver.dataset import NoCommunicationDataset


class Stage1(Stage):
    def _run(self):
        print('Doing some things')
        return 0


class Stage2(Stage1):
    previous_stages = Stage1


class Stage3(Stage1):
    previous_stages = (Stage1, Stage2)


def example_pipeline(d):
    final_stage = Stage3(dataset=d)
    final_stage.exit_status = 9
    luigi.build([final_stage], local_scheduler=True)
    return final_stage.exit_status


if __name__ == '__main__':
    cfg.load_config_file(join(dirname(dirname(__file__)), 'etc', 'example_analysisdriver.yaml'))
    dataset = NoCommunicationDataset('test_dataset', {'proc_id': 'test_dataset_15_02_2016_12:14:31'})
    print(example_pipeline(dataset))
