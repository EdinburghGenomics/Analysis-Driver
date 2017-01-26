import luigi
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
    previous_stages = ()
    exit_status = None

    dataset = EGCGParameter()

    @property
    def stage_name(self):
        if self.__stagename__:
            return self.__stagename__
        return self.__class__.__name__.lower()

    def requires(self):
        """
        Generates prior Stages from self.previous_stages, which should be either a single luigi.Task or a
        tuple. If it's a tuple, each element can be a luigi.Task or a
        tuple[luigi.Task,dict[str,luigi.Parameter]].
        """
        p = self.previous_stages
        if type(p) is not tuple:
            p = (p,)

        for s in p:
            if isinstance(s, type):
                s = s(dataset=self.dataset)
            yield s


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

    def _run(self):
        raise NotImplementedError


class RestAPITarget(luigi.Target):
    def __init__(self, stage):
        self.stage = stage

    @property
    def rest_api_stage(self):
        return self.stage.dataset.get_stage(self.stage.stage_name)

    def exists(self):
        s = self.rest_api_stage
        return s and bool(s.get('date_finished')) and s.get('exit_status') == 0


# example Luigi workflow below
from os.path import join, dirname
from time import sleep
from analysis_driver.dataset import NoCommunicationDataset


class ExampleStage(BasicStage):
    def run(self):
        print(self.stage_name)
        print('starting')
        sleep(5)
        print('done')
        return 0


class Stage1(ExampleStage):
    pass


class Stage2(ExampleStage):
    previous_stages = Stage1


class Stage3(ExampleStage):
    previous_stages = Stage1


class Stage4(ExampleStage):
    previous_stages = (Stage2, Stage3)


def example_pipeline(d):
    final_stage = Stage4(dataset=d)
    luigi.build([final_stage], local_scheduler=True)
    return final_stage.exit_status


if __name__ == '__main__':
    cfg.load_config_file(join(dirname(dirname(__file__)), 'etc', 'example_analysisdriver.yaml'))
    print(
        example_pipeline(
            NoCommunicationDataset(
                'test_dataset',
                {'proc_id': 'test_Dataset_15_02_2016_12:14:31'}
            )
        )
    )
