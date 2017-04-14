import os
import shutil
from luigi import Parameter, ListParameter
from egcg_core import executor
from analysis_driver.config import default as cfg
from analysis_driver.segmentation import Stage


class ContaminationCheck(Stage):
    fastq_files = ListParameter()

    def _fastqscreen_command(self):
        assert 1 <= len(self.fastq_files) <= 2, 'Bad number of fastqs: %s' % self.fastq_files
        return '%s --aligner bowtie2 %s --conf %s --force' % (
            cfg['tools']['fastqscreen'],
            ' '.join(self.fastq_files),
            cfg['contamination-check']['fastqscreen_conf']
        )

    @property
    def fastqscreen_expected_outfiles(self):
        return [f.rstrip('.fastq.gz') + '_screen.txt' for f in self.fastq_files]

    def _run(self):
        fastqscreen_run_command = self._fastqscreen_command()
        return executor.execute(
            fastqscreen_run_command,
            job_name='fastqscreen',
            working_dir=self.job_dir,
            cpus=2,
            mem=10
        ).join()


class VerifyBamID(Stage):
    bam_file = Parameter()

    @property
    def filtered_bam(self):
        return os.path.join(self.job_dir, self.dataset.name + '_chr22.bam')

    def _filter_bam(self):
        # use only chromosome 22 for speed
        return executor.execute(
            cfg['tools']['samtools'] + ' view -b %s chr22 > %s' % (self.bam_file, self.filtered_bam),
            job_name='filter_bam22',
            working_dir=self.job_dir,
            cpus=1,
            mem=2,
            log_commands=False
        ).join()

    def _index_filtered_bam(self):
        return executor.execute(
            cfg['tools']['samtools'] + ' index %s' % self.filtered_bam,
            job_name='index_bam22',
            working_dir=self.job_dir,
            cpus=1,
            mem=2
        ).join()

    def _verify_bam_id(self):
        cmd = '%s --bam %s --vcf %s --out %s' % (
            cfg['tools']['verifybamid'],
            self.filtered_bam,
            cfg['contamination-check']['population_vcf'],
            os.path.join(self.job_dir, self.dataset.name + '-chr22-vbi')
        )
        exit_status = executor.execute(
            cmd,
            job_name='verify_bam_id',
            working_dir=self.job_dir,
            cpus=1,
            mem=4
        ).join()
        sample_vbi_self = os.path.join(self.job_dir, self.dataset.name + '-chr22-vbi.selfSM')
        if os.path.exists(sample_vbi_self):
            bam_dir = os.path.dirname(self.bam_file)
            dest = os.path.join(bam_dir, os.path.basename(sample_vbi_self))
            shutil.copyfile(sample_vbi_self, dest)
        return exit_status

    def _run(self):
        exit_status = 0
        exit_status += self._filter_bam()
        exit_status += self._index_filtered_bam()
        exit_status += self._verify_bam_id()
        return exit_status


class VCFStats(Stage):
    vcf_file = Parameter()

    def _run(self):
        name, ext = os.path.splitext(self.vcf_file)
        stats_file = name + '.stats'
        cmd = '%s vcfstats %s > %s' % (cfg['tools']['rtg'], self.vcf_file, stats_file)
        exit_status = executor.execute(
            cmd,
            job_name='rtg_vcfstats',
            working_dir=self.job_dir,
            cpus=4,
            mem=32
        ).join()
        return exit_status
