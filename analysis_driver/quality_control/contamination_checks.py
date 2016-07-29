import os

import shutil

from analysis_driver import executor
from analysis_driver.config import default as cfg
from .quality_control_base import QualityControl


class ContaminationCheck(QualityControl):
    def __init__(self, dataset, working_dir, fastq_files):
        super().__init__(dataset, working_dir)
        self.fastq_files = fastq_files
        self.contamination_cfg = cfg.get('contamination-check')
        self.tools = cfg.get('tools')
        self.fastqscreen_expected_outfiles = None

    def _fastqscreen_command(self):
        """
        :return string: the command used to run fastq screen
        """
        if len(self.fastq_files) == 1:

            fastqscreen_bin = (self.tools.get('fastqscreen'))
            confPath = (self.contamination_cfg.get('fastqscreen_conf'))
            fastq1 = self.fastq_files[0]

            fastqscreen_command = ("{} --aligner bowtie2 {} --conf {} --force".format(fastqscreen_bin,
                                                                               fastq1,
                                                                               confPath))
            return fastqscreen_command

        elif len(self.fastq_files) == 2:
            fastqscreen_bin = (self.tools.get('fastqscreen'))
            confPath = (self.contamination_cfg.get('fastqscreen_conf'))
            fastq1 = self.fastq_files[0]
            fastq2 = self.fastq_files[1]

            fastqscreen_command = ("{} --aligner bowtie2 {} {} --conf {} --force".format(fastqscreen_bin,
                                                                                fastq1,
                                                                                fastq2,
                                                                                confPath))
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
        self.dataset.start_stage('run_fastqscreen')
        fastqscreen_executor = executor.execute(
            fastqscreen_run_command,
            job_name='fastqscreen',
            working_dir=self.working_dir,
            cpus=2,
            mem=10
        )
        exit_status = fastqscreen_executor.join()
        self.dataset.end_stage('run_fastqscreen', exit_status)
        return fastqscreen_expected_outfiles

    def run(self):
        try:
            self.fastqscreen_expected_outfiles = self._run_fastqscreen()
        except Exception as e:
            self.exception = e
            self.exit_status = 8

    def join(self, timeout=None):
        super().join(timeout=timeout)
        if self.exception:
            raise self.exception
        return self.fastqscreen_expected_outfiles


class VerifyBamId(QualityControl):
    """This class runs verifyBamId on a subset of the provided bam file"""
    def __init__(self, dataset, working_dir, bam_file):
        super().__init__(dataset, working_dir)
        self.input_bam = bam_file
        self.exit_status = None

    def _filter_bam(self):
        self.filtered_bam = os.path.join(self.working_dir, self.dataset.name + '_chr22.bam')
        cmd = cfg.query('tools', 'samtools') + ' view -b %s chr22 > %s' % (self.input_bam, self.filtered_bam)
        return executor.execute(
                cmd,
                job_name='filter_bam22',
                working_dir=self.working_dir,
                cpus=1,
                mem=2,
                log_commands=False
        ).join()

    def _index_filtered_bam(self):
        cmd = cfg.query('tools', 'samtools') + ' index %s' % (self.filtered_bam)
        return executor.execute(
                cmd,
                job_name='index_bam22',
                working_dir=self.working_dir,
                cpus=1,
                mem=2
        ).join()

    def _verify_bam_id(self):
        population_vcf = cfg.query('contamination-check', 'population_vcf')
        cmd = '%s --bam %s --vcf %s --out %s' % (
            cfg.query('tools', 'verifybamid'),
            self.filtered_bam,
            population_vcf,
            os.path.join(self.working_dir, self.dataset.name + '-chr22-vbi')
        )
        exit_status = executor.execute(
            cmd,
            job_name='verify_bam_id',
            working_dir=self.working_dir,
            cpus=1,
            mem=4
        ).join()
        sample_vbi_self = os.path.join(self.working_dir, self.dataset.name + '-chr22-vbi.selfSM')
        if os.path.exists(sample_vbi_self):
            bam_dir = os.path.dirname(self.input_bam)
            dest = os.path.join(bam_dir, os.path.basename(sample_vbi_self))
            shutil.copyfile(sample_vbi_self, dest)
        return exit_status

    def _contamination_check(self):
        exit_status = 0
        self.dataset.start_stage('verify_bam_id')
        exit_status += self._filter_bam()
        exit_status += self._index_filtered_bam()
        exit_status += self._verify_bam_id()
        self.dataset.end_stage('verify_bam_id', exit_status)
        return exit_status

    def run(self):
        try:
            self.exit_status = self._contamination_check()
        except Exception as e:
            self.exception = e

    def join(self, timeout=None):
        super().join(timeout=timeout)
        if self.exception:
            raise self.exception
        return self.exit_status


class VCFStats(QualityControl):
    """This class runs vcfstats from rtg on a filtered vcf file"""
    def __init__(self, dataset, working_dir, vcf_file):
        super().__init__(dataset, working_dir)
        self.vcf_file = vcf_file
        self.exit_status = None

    def _vcf_stats(self):
        name, ext = os.path.splitext(self.vcf_file)
        stats_file = name + '.stats'
        cmd = '%s vcfstats %s > %s' % (cfg.query('tools', 'rtg'), self.vcf_file, stats_file)
        exit_status = executor.execute(
            cmd,
            job_name='rtg_vcfstats',
            working_dir=self.working_dir,
            cpus=4,
            mem=32
        ).join()
        return exit_status

    def run(self):
        try:
            self.exit_status = self._vcf_stats()
        except Exception as e:
            self.exception = e

    def join(self, timeout=None):
        super().join(timeout=timeout)
        if self.exception:
            raise self.exception
        return self.exit_status


