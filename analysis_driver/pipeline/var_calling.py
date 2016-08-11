from os.path import join
from analysis_driver import quality_control as qc
from analysis_driver.pipeline import Stage, gatk
from analysis_driver.pipeline.common import MergeFastqs, Fastqc, CoverageStats
from analysis_driver.config import default as cfg
from analysis_driver.exceptions import PipelineError
from analysis_driver.util import bash_commands
from egcg_core import executor, util


class VarCallingStage(Stage):
    @property
    def fastq_pair(self):
        return util.find_files(self.job_dir, 'merged', self.user_sample_id + '_R?.fastq.gz')

    @property
    def output_bam(self):
        return join(self.job_dir, self.dataset_name + '.bam')

    def _run(self):
        raise NotImplementedError


class SpeciesContaminationCheck(VarCallingStage):
    previous_stages = MergeFastqs

    def _run(self):
        species_contamination_check = qc.ContaminationCheck(self.dataset, self.job_dir, [self.fastq_pair[0]])
        species_contamination_check.start()
        species_contamination_check.join()
        return species_contamination_check.exit_status


class BWAAlignment(VarCallingStage):
    previous_stages = MergeFastqs

    def _run(self):
        reference = cfg.query('references', self.dataset.species, 'fasta')
        if not reference:
            raise PipelineError('Could not find reference for species %s in sample %s ' % (self.dataset.species, self.dataset_name))

        return executor.execute(
            bash_commands.bwa_mem_samblaster(
                self.fastq_pair,
                reference,
                self.output_bam,
                read_group={'ID': '1', 'SM': self.dataset.user_sample_id, 'PL': 'illumina'},
                thread=16
            ),
            job_name='bwa_mem',
            working_dir=self.job_dir,
            cpus=16,
            mem=64
        ).join()


class BamtoolsStats(VarCallingStage):
    previous_stages = BWAAlignment

    def _run(self):
        return executor.execute(
            bash_commands.bamtools_stats(self.output_bam, join(self.job_dir, 'bamtools_stats.txt')),
            job_name='bamtools',
            working_dir=self.job_dir,
            cpus=1,
            mem=4,
            log_commands=False
        ).join()


class BasicQC(VarCallingStage):
    @property
    def previous_stages(self):
        return Fastqc(previous_stages=MergeFastqs, fastqs=self.fastq_pair), SpeciesContaminationCheck, BamtoolsStats, CoverageStats(bam_file=self.output_bam)

    def _run(self):
        return 0


class VarCalling(VarCallingStage):
    previous_stages = (BasicQC, gatk.Tabix)

    def _run(self):
        return 0
