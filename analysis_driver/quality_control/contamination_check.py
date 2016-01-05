__author__ = 'mwham'
from threading import Thread
import os.path
from analysis_driver.notification import default as ntf
from analysis_driver.writer import bash_commands
from analysis_driver import executor
from analysis_driver.app_logging import AppLogger


class ContaminationCheck(AppLogger, Thread):
    def __init__(self, run_id, file_prefix):
        Thread.__init__(self)
        self.run_id = run_id
        self.bam_file = file_prefix + '.bam'
        self.vcf_file = file_prefix + '.vcf.gz'  # TODO: clean up input/output file locations
        self.exception = None

    def _contamination_check(self):
        exit_status = 0
        exit_status += self._verify_bam_id()
        return exit_status

    def _verify_bam_id(self):
        ntf.start_stage('verify_bam_id')
        cmd = bash_commands.verify_bam_id(
            self.bam_file,
            self.vcf_file,
            out_prefix=os.path.join(os.path.dirname(self.vcf_file), 'verify_bam_id_' + self.run_id)
        )
        exit_status = executor.execute(
            [cmd],
            job_name='verify_bam_id',
            run_id=self.run_id,
            cpus=4,
            mem=8
        ).join()
        ntf.end_stage('verify_bam_id', exit_status)

    def run(self):
        try:
            self._contamination_check()
        except Exception as e:
            self.exception = e

    def join(self, timeout=None):
        super().join(timeout=timeout)
        if self.exception:
            raise self.exception
