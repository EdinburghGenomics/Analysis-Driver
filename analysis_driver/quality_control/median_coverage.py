import os
from threading import Thread
from analysis_driver import executor
from analysis_driver.app_logging import AppLogger
from analysis_driver.config import default as cfg
from analysis_driver.notification import default as ntf

class GatkDepthofCoverage(AppLogger, Thread):
    def __init__(self, sample_id, bam_file):
        self.sample_id = sample_id
        self.bam_file = bam_file
        self.work_directory = os.path.join(cfg['jobs_dir'], self.sample_id)
        self.exception = None
        Thread.__init__(self)

    def _get_gatk_depthofcoverage_command(self):
        reference = cfg['bcbio'] + 'genomes/Hsapiens/' + cfg['genome'] + '/seq/' + cfg['genome'] + '.fa'
        gatk_depthofcoverage_out_file = ((self.bam_file.rstrip('bam')) + 'depthofcoverage')
        gatk_depthofcoverage_command = 'java -jar GenomeAnalysisTK.jar' \
                                       ' -T DepthOfCoverage ' \
                                       '-R %s ' \
                                       '-o %s ' \
                                       '-I %s' % (reference, gatk_depthofcoverage_out_file, self.bam_file)
        return gatk_depthofcoverage_command, gatk_depthofcoverage_out_file

    def _run_gatk_depthofcoverage(self):
        """
        :return string: the expected outfile from GATK depthofcoverage
        """
        [depthofcoverage_command], depthofcoverage_out_file = self._get_gatk_depthofcoverage_command()
        ntf.start_stage('run_gatk_depthofcoverage')
        depthofcoverage_executor = executor.execute(
            depthofcoverage_command,
            job_name='depthofcoverage',
            working_dir=self.work_directory,
            cpus=2,
            mem=10
        )
        exit_status = depthofcoverage_executor.join()
        ntf.end_stage('run_gatk_depthofcoverage', exit_status)
        return depthofcoverage_out_file

    def run(self):
        try:
            self._run_gatk_depthofcoverage()
        except Exception as e:
            self.exception = e

    def join(self, timeout=None):
        super().join(timeout=timeout)
        if self.exception:
            raise self.exception
        return self.median_coverage_expected_outfiles