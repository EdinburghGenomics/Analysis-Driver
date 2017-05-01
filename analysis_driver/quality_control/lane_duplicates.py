import os
from luigi import Parameter
from egcg_core import executor
from analysis_driver.config import default as cfg
from analysis_driver.segmentation import Stage

class WellDuplicates(Stage):
    run_directory = Parameter()
    output_directory = Parameter()

    def _welldups_cmd(self):
        output_file = os.path.join(self.output_directory, self.dataset.name + '.wellduplicate')
        return '{welldups} -f {coords} -r {run_dir} -s hiseq_x > {outfile} 2> {outfile}.err'.format(
            welldups=cfg.query('tools', 'well_duplicate'),
            coords=cfg['well_duplicate']['coord_file'], run_dir=self.run_directory, outfile=output_file
        )

    def _run(self):
        return executor.execute(
            self._welldups_cmd(),
            job_name='welldup', working_dir=self.job_dir, cpus=1, mem=2, log_commands=False
        ).join()

