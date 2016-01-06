__author__ = 'mwham'
from threading import Thread
import os.path
from analysis_driver.notification import default as ntf
from analysis_driver.writer import bash_commands
from analysis_driver import executor
from analysis_driver.app_logging import AppLogger


class ContaminationCheck(AppLogger, Thread):
    def __init__(self, run_id, sample_id, input_dir, output_dir):
        Thread.__init__(self)
        self.run_id = run_id
        self.sample_id = sample_id

        self.input_bam = os.path.join(input_dir, sample_id + '.bam')
        self.input_vcf = os.path.join(input_dir, sample_id + '.vcf.gz')
        self.automsomal_vcf = None
        self.output_dir = output_dir
        self.exception = None

    def _contamination_check(self):
        exit_status = 0
        exit_status += self._verify_bam_id()
        return exit_status

    def _remove_non_autosomes(self):
        self.automsomal_vcf = os.path.join(self.output_dir, self.sample_id + '_autosomes.vcf.gz')
        cmd = bash_commands.remove_non_autosomes(self.input_vcf, self.automsomal_vcf)
        executor.execute([cmd], job_name='remove_non_autosomes', run_id=self.run_id, cpus=1, mem=2).join()

    def _verify_bam_id(self):
        ntf.start_stage('verify_bam_id')
        if os.path.isfile(self.automsomal_vcf):
            cmd = bash_commands.verify_bam_id(
                self.input_bam,
                vcf_file=self.automsomal_vcf,
                out_prefix=self.output_dir + self.sample_id
            )
            exit_status = executor.execute(
                [cmd],
                job_name='verify_bam_id',
                run_id=self.run_id,
                cpus=4,
                mem=8
            ).join()
        else:
            self.error('No autosomal VCF found')
            exit_status = 1
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
