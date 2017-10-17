import os
from luigi import Parameter
from egcg_core import executor
from egcg_core.util import find_file
from luigi.parameter import ListParameter

from analysis_driver.tool_versioning import toolset
from analysis_driver.segmentation import Stage
from analysis_driver.util import bash_commands


class SamtoolsDepth(Stage):
    bam_file = Parameter()

    @property
    def samtools_depth_out_file(self):
        return os.path.splitext(find_file(self.bam_file))[0] + '.depth'

    def _run(self):
        self.info('/Generating depth file: %s', self.samtools_depth_out_file)
        return executor.execute(
            bash_commands.samtools_depth_command(
                self.job_dir,
                find_file(self.bam_file),
                self.samtools_depth_out_file
            ),
            job_name='samtoolsdepth',
            working_dir=self.job_dir,
            cpus=1,
            mem=6
        ).join()
