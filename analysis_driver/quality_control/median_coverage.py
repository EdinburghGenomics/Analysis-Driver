import os
from analysis_driver import executor
from analysis_driver.config import default as cfg
from analysis_driver.notification import default as ntf
from .quality_control_base import QualityControl


class SamtoolsDepth(QualityControl):
    def __init__(self, dataset, working_dir, bam_file):
        super().__init__(dataset, working_dir)
        self.bam_file = bam_file
        self.working_dir =  working_dir
        self.exception = None


    def _get_samtools_depth_command(self):
        samtools_bin = cfg['tools']['samtools']
        samtools_depth_out_file = ((self.bam_file).rstrip('bam') + 'depth')
        samtools_depth_command = "%s depth -a -a -q 0 -Q 0 %s > %s" % (samtools_bin, self.bam_file, samtools_depth_out_file)
        return samtools_depth_command, samtools_depth_out_file



    def _run_samtools_depth(self):
        """
        :return string: the expected outfile from samtools depth
        """
        samtools_depth_command, samtools_depth_out_file = self._get_samtools_depth_command()
        ntf.start_stage('run_samtools_depth')
        samtools_depth_executor = executor.execute(
            [samtools_depth_command],
            job_name='samtoolsdepth',
            working_dir=self.working_dir,
            cpus=2,
            mem=10
        )
        exit_status = samtools_depth_executor.join()
        ntf.end_stage('run_samtools_depth', exit_status)
        return samtools_depth_out_file

    def _get_depth_histogram_command(self):
        samtools_depth_outfile = self._run_samtools_depth()
        depth_histogram_outfile = samtools_depth_outfile.rstrip('depth') + 'hist'
        depth_histogram_command = "awk -F '\t' '{print $3}' %s | sort | uniq -c | sort -nr | sort -k2 -n > %s" % (samtools_depth_outfile,
                                                                                                    depth_histogram_outfile)
        return depth_histogram_command, depth_histogram_outfile

    def _run_depth_histogram_command(self):
        depth_histogram_command, depth_histogram_outfile = self._get_depth_histogram_command()
        ntf.start_stage('run_depth_histogram')
        depth_histogram_executor = executor.execute(
            [depth_histogram_command],
            job_name='depthhistogram',
            working_dir=self.working_dir,
            cpus=2,
            mem=10
        )
        exit_status = depth_histogram_executor.join()
        ntf.end_stage('run_depth_histogram', exit_status)
        return depth_histogram_outfile



    def run(self):
        try:
            self.median_coverage_expected_outfiles = self._run_depth_histogram_command()
        except Exception as e:
            self.exception = e

    def join(self, timeout=None):
        super().join(timeout=timeout)
        if self.exception:
            raise self.exception
        return self.median_coverage_expected_outfiles