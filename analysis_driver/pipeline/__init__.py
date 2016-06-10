import luigi
from os.path import join
from analysis_driver.config import default as cfg
from analysis_driver.exceptions import PipelineError


class Stage(luigi.Task):
    __stagename__ = None
    exit_status = None
    previous_stages = ()
    dataset = luigi.Parameter()

    @property
    def dataset_name(self):
        return self.dataset.name

    @property
    def job_dir(self):
        return join(cfg['jobs_dir'], self.dataset_name)

    @property
    def input_dir(self):
        return join(cfg.get('intermediate_dir', cfg['input_dir']), self.dataset_name)

    def output(self):  # if <stage_name>.done is present, the stage is complete
        return luigi.LocalTarget(self._stage_lock_file())

    def run(self):
        self.dataset.start_stage(self.stage_name)
        self.exit_status = self._run()
        self.dataset.end_stage(self.stage_name, self.exit_status)
        if self.exit_status == 0:
            self._touch(self._stage_lock_file())
        else:
            raise PipelineError('Exit status was %s. Stopping' % self.exit_status)

    @property
    def stage_name(self):
        if self.__stagename__:
            return self.__stagename__
        return self.__class__.__name__.lower()

    def _run(self):
        raise NotImplementedError

    @staticmethod
    def _touch(f):
        open(f, 'w').close()

    def _stage_lock_file(self):
        return join(self.job_dir, '.' + self.stage_name + '.done')

    def requires(self):
        if len(self.previous_stages) == 1:
            return self.previous_stages[0](dataset=self.dataset)
        else:
            return [s(dataset=self.dataset) for s in self.previous_stages]
