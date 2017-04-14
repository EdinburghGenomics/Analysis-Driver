import os
from luigi import Parameter
from egcg_core import executor
from analysis_driver.config import default as cfg
from analysis_driver.segmentation import Stage


class WellDuplicates(Stage):
    run_directory = Parameter()
    output_directory = Parameter()

    def _run(self):
        output_file = os.path.join(self.output_directory, self.dataset.name + '.wellduplicate')
        output_err = output_file + '.err'
        coord_file = cfg.query('well_duplicate', 'coord_file')

        cmd = cfg.query('tools', 'well_duplicate') + ' -f %s -r %s -s hiseq_x > %s 2> %s' % (
            coord_file, self.run_directory, output_file, output_err
        )
        return executor.execute(
            cmd, job_name='welldup', working_dir=self.job_dir, cpus=1, mem=2, log_commands=False
        ).join()
