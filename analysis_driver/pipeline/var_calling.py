from os.path import join
from analysis_driver import quality_control as qc
from analysis_driver.pipeline import Stage
from analysis_driver.pipeline import gatk
from analysis_driver.transfer_data import prepare_sample_data
from analysis_driver.driver import _bcbio_prepare_sample
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


class MergeFastqs(VarCallingStage):
    def _run(self):
        fastq_files = prepare_sample_data(self.dataset)
        _bcbio_prepare_sample(self.job_dir, self.dataset_name, fastq_files)
        self.debug('Sample fastq files: ' + str(self.fastq_pair))


class MergedFastqc(VarCallingStage):
    previous_stages = MergeFastqs

    def _run(self):
        return executor.execute(
            *[bash_commands.fastqc(fastq_file) for fastq_file in self.fastq_pair],
            job_name='fastqc2',
            working_dir=self.job_dir,
            cpus=1,
            mem=2
        ).join()


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


class CoverageStats(VarCallingStage):
    previous_stages = BWAAlignment

    def _run(self):
        coverage_statistics_histogram = qc.SamtoolsDepth(self.dataset, self.job_dir, self.output_bam)
        coverage_statistics_histogram.start()
        coverage_statistics_histogram.join()
        return coverage_statistics_histogram.exit_status


class BasicQC(VarCallingStage):
    previous_stages = (MergedFastqc, SpeciesContaminationCheck, BamtoolsStats, CoverageStats)

    def _run(self):
        return 0


class VarCalling(VarCallingStage):
    previous_stages = (MergedFastqc, SpeciesContaminationCheck, BamtoolsStats, CoverageStats, gatk.Tabix)

    def _run(self):
        return 0
