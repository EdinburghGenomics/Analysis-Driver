import os
from luigi import Parameter
from egcg_core import executor
from analysis_driver.config import default as cfg
from analysis_driver.segmentation import Stage
from analysis_driver.tool_versioning import toolset


class WellDuplicates(Stage):
    run_directory = Parameter()
    output_directory = Parameter()

    def _welldups_cmd(self):
        output_file = os.path.join(self.output_directory, self.dataset.name + '.wellduplicate')
        return '{welldups} -f {coords} -r {run_dir} -s hiseq_x > {outfile} 2> {outfile}.err'.format(
            welldups=toolset['well_duplicates'],
            coords=cfg['well_duplicates']['coord_file'], run_dir=self.run_directory, outfile=output_file
        )

    def _run(self):
        return executor.execute(
            self._welldups_cmd(),
            job_name='welldup', working_dir=self.job_dir, cpus=1, mem=2, log_commands=False
        ).join()
