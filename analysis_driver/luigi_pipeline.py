import luigi
from os.path import join
from analysis_driver.exceptions import PipelineError
from analysis_driver.app_logging import logging_default as log_cfg
from analysis_driver.config import default as cfg


app_logger = log_cfg.get_logger('luigi')


class Stage(luigi.Task):
    dataset = luigi.Parameter()
    
    @property
    def dataset_name(self):
        return self.dataset.name

    @property
    def job_dir(self):
        return join(cfg['jobs_dir'], self.dataset_name)

    def output(self):  # if <stage_name>.done is present, the stage is complete
        return luigi.LocalTarget(self._stage_lock_file())

    def run(self):
        self.dataset.start_stage(self.stage_name)
        exit_status = self._run()
        self.dataset.end_stage(self.stage_name, exit_status)
        if exit_status == 0:
            self._touch(self._stage_lock_file())
        else:
            raise PipelineError('Exit status was %s. Stopping' % exit_status)

    @property
    def stage_name(self):
        return self.__class__.__name__.lower()

    def _run(self):
        raise NotImplementedError

    @staticmethod
    def _touch(f):
        open(f, 'w').close()

    def _stage_lock_file(self):
        return join(self.job_dir, '.' + self.stage_name + '.done')


def pipeline(dataset):

    luigi.run(
        cmdline_args=[
            '--dataset', dataset,
        ],
        # main_task_cls=DataOutput,
        local_scheduler=True
    )
