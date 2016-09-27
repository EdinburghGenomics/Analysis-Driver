import os
from egcg_core import executor
from analysis_driver.config import default as cfg
from analysis_driver.notification import default as ntf
from .quality_control_base import QualityControl


class SamtoolsDepth(QualityControl):
    def __init__(self, dataset, working_dir, bam_file):
        super().__init__(dataset, working_dir)
        self.bam_file = bam_file
        self.working_dir = working_dir

    def _get_samtools_depth_command(self):
        samtools_bin = cfg['tools']['samtools']
        name, ext = os.path.splitext(self.bam_file)
        samtools_depth_out_file = name + '.depth'
        samtools_depth_command = """%s depth -a -a -q 0 -Q 0 %s | awk -F "\t" '{array[$1"\t"$3]+=1} END{for (val in array){print val"\t"array[val]}}' | sort -k 1,1 -nk 2,2 > %s""" % (samtools_bin, self.bam_file, samtools_depth_out_file)
        return samtools_depth_command, samtools_depth_out_file

    def _run_samtools_depth(self):
        """
        :return string: the expected outfile from samtools depth
        """
        samtools_depth_command, samtools_depth_out_file = self._get_samtools_depth_command()
        ntf.start_stage('run_samtools_depth')
        samtools_depth_executor = executor.execute(
            samtools_depth_command,
            job_name='samtoolsdepth',
            working_dir=self.working_dir,
            cpus=1,
            mem=6
        )
        exit_status = samtools_depth_executor.join()
        ntf.end_stage('run_samtools_depth', exit_status)
        return samtools_depth_out_file


    def run(self):
        try:
            self.median_coverage_expected_outfiles = self._run_samtools_depth()
        except Exception as e:
            self.exception = e

    def join(self, timeout=None):
        super().join(timeout=timeout)
        if self.exception:
            raise self.exception
        return self.median_coverage_expected_outfiles