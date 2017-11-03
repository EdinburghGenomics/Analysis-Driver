import os
from luigi import Parameter
from egcg_core import executor
from egcg_core.util import find_file
from analysis_driver.tool_versioning import toolset
from analysis_driver.segmentation import Stage


class SamtoolsDepth(Stage):
    bam_file = Parameter()

    @property
    def samtools_depth_out_file(self):
        return os.path.splitext(find_file(self.bam_file))[0] + '.depth'

    def _samtools_depth_command(self):
        return (
            '%s depth -a -a -q 0 -Q 0 %s | '
            'awk -F "\t" \'{array[$1"\t"$3]+=1} END{for (val in array){print val"\t"array[val]}}\' | '
            'sort -T %s -k 1,1 -nk 2,2 > %s'
        ) % (toolset['samtools'], find_file(self.bam_file), self.job_dir, self.samtools_depth_out_file)

    def _run(self):
        self.info('/Generating depth file: %s', self.samtools_depth_out_file)
        return executor.execute(
            self._samtools_depth_command(),
            job_name='samtoolsdepth',
            working_dir=self.job_dir,
            cpus=1,
            mem=6
        ).join()
