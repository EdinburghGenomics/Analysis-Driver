import os
from threading import Thread
from analysis_driver.config import default as cfg
from analysis_driver import executor
from analysis_driver.notification import default as ntf
from analysis_driver.app_logging import AppLogger



class ContaminationCheck(AppLogger, Thread):

    def __init__(self, fastq_files, working_dir):
        self.fastq_files = fastq_files
        self.working_dir = working_dir
        self.contamination_cfg = cfg.get('contamination-check')
        self.tools = cfg.get('tools')
        Thread.__init__(self)
        self.exception = None

    def _fastqscreen_command(self):
        """
        :return list: the command used to run fastq screen
        """
        if len(self.fastq_files) == 1:

            fastqscreen_bin = (self.tools.get('fastqscreen'))
            confPath = (self.contamination_cfg.get('fastqscreen_conf'))
            fastq1 = self.fastq_files[0]

            fastqscreen_command = ["{} --aligner bowtie2 {} --conf {} --force".format(fastqscreen_bin,
                                                                               fastq1,
                                                                               confPath)]
            return fastqscreen_command

        elif len(self.fastq_files) == 2:
            fastqscreen_bin = (self.tools.get('fastqscreen'))
            confPath = (self.contamination_cfg.get('fastqscreen_conf'))
            fastq1 = self.fastq_files[0]
            fastq2 = self.fastq_files[1]

            fastqscreen_command = ["{} --aligner bowtie2 {} {} --conf {} --force".format(fastqscreen_bin,
                                                                                fastq1,
                                                                                fastq2,
                                                                                confPath)]
            return fastqscreen_command
        else:
           raise ValueError('Bad number of fastqs: ' + str(self.fastq_files))

    def _get_expected_outfiles(self):
        if len(self.fastq_files) == 1:
            expected_outfiles = [((self.fastq_files[0]).rstrip('.fastq.gz') + '_screen.txt')]
            return expected_outfiles
        elif len(self.fastq_files) == 2:
            expected_outfiles = [((self.fastq_files[0]).rstrip('.fastq.gz') + '_screen.txt'), ((self.fastq_files[1]).rstrip('.fastq.gz') + '_screen.txt')]
            return expected_outfiles
        else:
            raise ValueError('Bad number of fastqs: ' + str(self.fastq_files))


    def _run_fastqscreen(self):
        """
        :return list: a list of the expected outfiles from fastq screen
        """
        fastqscreen_run_command = self._fastqscreen_command()
        fastqscreen_expected_outfiles = self._get_expected_outfiles()
        ntf.start_stage('run_fastqscreen')
        fastqscreen_executor = executor.execute(
            fastqscreen_run_command,
            job_name='fastqscreen',
            working_dir=self.working_dir,
            cpus=2,
            mem=10
        )
        exit_status = fastqscreen_executor.join()
        ntf.end_stage('run_fastqscreen', exit_status)
        return fastqscreen_expected_outfiles

    def run(self):
        try:
            self.fastqscreen_expected_outfiles = self._run_fastqscreen()
        except Exception as e:
            self.exception = e

    def join(self, timeout=None):
        super().join(timeout=timeout)
        if self.exception:
            raise self.exception
        return self.fastqscreen_expected_outfiles

