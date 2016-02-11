import os
from analysis_driver.config import default as cfg
from analysis_driver import executor
from analysis_driver.notification import default as ntf


class ContaminationCheck():

    def __init__(self, fastq_files, sample_id):
        self.fastq_files = fastq_files
        self.sample_id = sample_id
        self.contamination_cfg = cfg.get('contamination-check')

    def _fastqscreen_command(self):
        """
        :return list: the command used to run fastq screen
        """
        if len(self.fastq_files) == 1:

            fastqscreen_bin = (self.contamination_cfg.get('fastq_screen_bin').rstrip('/') + '/')
            confPath = (self.contamination_cfg.get('fastq_screen_conf').rstrip('/') + '/')
            fastq1 = self.fastq_files[0]

            fastqscreen_command = ["{} --aligner bowtie2 {} --conf {}".format(fastqscreen_bin,
                                                                               fastq1,
                                                                               confPath)]
            return fastqscreen_command

        elif len(self.fastq_files) == 2:
            fastqscreen_bin = (self.contamination_cfg.get('fastq_screen_bin').rstrip('/') + '/')
            confPath = (self.contamination_cfg.get('fastq_screen_conf').rstrip('/') + '/')
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
        fastqscreen_expected_outfiles = self._get_expected_outfiles()
        ntf.start_stage('fastqscreen_contamination_check')
        kontaminant_executor = executor.execute(
            fastqscreen_run_command,
            job_name='fastqscreen',
            sample_id=self.sample_id,
            cpus=4,
            mem=25
        )
        exit_status = kontaminant_executor.join()
        ntf.end_stage('fastqscreen_contamination_check', exit_status)
        return fastqscreen_expected_outfiles
