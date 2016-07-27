import os
from egcg_core import executor
from .quality_control_base import QualityControl
from analysis_driver.config import default as cfg


class WellDuplicates(QualityControl):
    """

    """
    def __init__(self, dataset, working_dir, output_directory, run_directory ):
        super().__init__(dataset, working_dir)
        self.run_directory = run_directory
        self.output_directory = output_directory
        self.working_dir = working_dir

    def _well_duplicates(self):
        """
        """
        output_file = os.path.join(self.output_directory, self.dataset.name + '.wellduplicate')
        output_err = os.path.join(self.output_directory, self.dataset.name + '.wellduplicate.err')
        coord_file = cfg.query('well_duplicate', 'coord_file')

        cmd = cfg.query('tools', 'well_duplicate') + " -f %s -r %s -s hiseq_x > %s 2> %s"%(
            coord_file,
            self.run_directory,
            output_file,
            output_err
        )
        return executor.execute(
                cmd,
                job_name='welldup',
                working_dir=self.working_dir,
                cpus=1,
                mem=2,
                log_commands=False
        ).join()

    def run(self):
        try:
            self.dataset.start_stage('wellduplicate')
            self.exit_status = self._well_duplicates()
        except Exception as e:
            self.exception = e
            self.exit_status = 9
        finally:
            self.dataset.end_stage('wellduplicate', self.exit_status)

    def join(self, timeout=None):
        super().join(timeout=timeout)
        if self.exception:
            raise self.exception
        return self.exit_status
