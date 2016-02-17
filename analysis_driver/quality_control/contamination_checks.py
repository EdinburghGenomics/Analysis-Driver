import os
from threading import Thread
from analysis_driver.config import default as cfg
from analysis_driver import executor
from analysis_driver.notification import default as ntf
from analysis_driver.app_logging import AppLogger



class ContaminationCheck(AppLogger, Thread):

    def __init__(self, fastq_files, sample_id):
        self.fastq_files = fastq_files
        self.sample_id = sample_id
        self.contamination_cfg = cfg.get('contamination-check')
        Thread.__init__(self)
        self.exception = None

    def _fastqscreen_command(self):
        """
        :return list: the command used to run fastq screen
        """
        if len(self.fastq_files) == 1:

            fastqscreen_bin = (self.contamination_cfg.get('fastq_screen_bin'))
            confPath = (self.contamination_cfg.get('fastq_screen_conf'))
            fastq1 = self.fastq_files[0]

            fastqscreen_command = ["{} --aligner bowtie2 {} --conf {}".format(fastqscreen_bin,
                                                                               fastq1,
                                                                               confPath)]
            return fastqscreen_command

        elif len(self.fastq_files) == 2:
            fastqscreen_bin = (self.contamination_cfg.get('fastq_screen_bin'))
            confPath = (self.contamination_cfg.get('fastq_screen_conf'))
            fastq1 = self.fastq_files[0]
            fastq2 = self.fastq_files[1]

            fastqscreen_command = ["{} --aligner bowtie2 {} {} --conf {}".format(fastqscreen_bin,
                                                                                fastq1,
                                                                                fastq2,
                                                                                confPath)]
            return fastqscreen_command
        else:
           raise ValueError('Bad number of fastqs: ' + str(self.fastq_files))

    def _get_expected_outfiles(self):
        if len(self.fastq_files) == 1:
            expected_outfiles = [((self.fastq_files[0]).rstrip('.fastq') + '_screen.txt')]
            return expected_outfiles
        elif len(self.fastq_files) == 2:
            expected_outfiles = [((self.fastq_files[0]).rstrip('.fastq') + '_screen.txt'), ((self.fastq_files[1]).rstrip('.fastq') + '_screen.txt')]
            return expected_outfiles
        else:
            raise ValueError('Bad number of fastqs: ' + str(self.fastq_files))

    def _run_fastqscreen(self):
        fastqscreen_run_command = self._fastqscreen_command()
        print('printing fastqscreen run command')
        print(fastqscreen_run_command)
        print('\n')
        fastqscreen_expected_outfiles = self._get_expected_outfiles()
        ntf.start_stage('fastqscreen_contamination_check')
        fastqscreen_executor = executor.execute(
            fastqscreen_run_command,
            job_name='fastqscreen',
            run_id=self.sample_id,
            cpus=2,
            mem=10
        )
        exit_status = fastqscreen_executor.join()
        ntf.end_stage('fastqscreen_contamination_check', exit_status)
        print('finished now')
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

