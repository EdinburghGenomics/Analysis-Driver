from analysis_driver import quality_control as qc
from analysis_driver.pipeline import Stage
from analysis_driver.transfer_data import prepare_sample_data
from analysis_driver.driver import _bcbio_prepare_sample, _run_bcbio, _link_results_files, _output_data
from analysis_driver.util import bash_commands
from egcg_core import executor, util


class BCBioStage(Stage):
    @property
    def fastq_pair(self):
        return util.find_files(self.job_dir, 'merged', self.user_sample_id + '_R?.fastq.gz')

    def _run(self):
        raise NotImplementedError


class MergeFastqs(BCBioStage):
    def _run(self):
        fastq_files = prepare_sample_data(self.dataset)
        _bcbio_prepare_sample(self.job_dir, self.dataset_name, fastq_files)
        self.debug('sample fastq files: ' + str(self.fastq_pair))


class MergedFastqc(BCBioStage):
    previous_stages = MergeFastqs

    def _run(self):
        return executor.execute(
            *[bash_commands.fastqc(fastq_file) for fastq_file in self.fastq_pair],
            job_name='fastqc2',
            working_dir=self.job_dir,
            cpus=1,
            mem=2
        ).join()


class GenotypeValidation(BCBioStage):
    previous_stages = MergeFastqs

    def _run(self):
        genotype_validation = qc.GenotypeValidation(self.dataset, self.job_dir, self.fastq_pair)
        genotype_validation.start()
        geno_valid_vcf_file, geno_valid_results = genotype_validation.join()
        self.info('Written files: ' + str(geno_valid_vcf_file) + ' ' + str(geno_valid_results))
        return genotype_validation.exit_status


class SpeciesContaminationCheck(BCBioStage):
    previous_stages = MergeFastqs

    def _run(self):
        species_contamination_check = qc.ContaminationCheck(self.dataset, self.job_dir, [self.fastq_pair[0]])
        species_contamination_check.start()
        species_contamination_check.join()
        return species_contamination_check.exit_status


class BCBio(BCBioStage):
    previous_stages = MergeFastqs

    def _run(self):
        return _run_bcbio(self.dataset_name, self.job_dir, self.fastq_pair).join()


class GenderValidation(BCBioStage):
    previous_stages = BCBio

    def _run(self):
        vcf_file = util.find_file(
            'samples_%s-merged' % self.dataset_name,
            'final',
            '*_' + self.dataset.user_sample_id,
            self.dataset.user_sample_id + '-joint-gatk-haplotype-joint.vcf.gz'
        )
        gender_validation = qc.GenderValidation(self.dataset, self.job_dir, vcf_file)
        gender_validation.start()
        return gender_validation.join()


class SampleContaminationCheck(BCBioStage):
    previous_stages = BCBio

    def _run(self):
        bam_file = util.find_file(
            'samples_%s-merged' % self.dataset_name,
            'final',
            self.dataset.user_sample_id,
            self.dataset.user_sample_id + '-ready.bam'
        )
        sample_contam = qc.VerifyBamId(self.dataset, self.job_dir, bam_file)
        sample_contam.start()
        return sample_contam.join()


class CoverageStats(BCBioStage):
    previous_stages = BCBio

    def _run(self):
        bam_file = util.find_file(
            'samples_%s-merged' % self.dataset_name,
            'final',
            self.dataset.user_sample_id,
            self.dataset.user_sample_id + '-ready.bam'
        )
        coverage_statistics_histogram = qc.SamtoolsDepth(self.dataset, self.job_dir, bam_file)
        coverage_statistics_histogram.start()
        coverage_statistics_histogram.join()
        return coverage_statistics_histogram.exit_status


class DataOutput(BCBioStage):
    previous_stages = (MergedFastqc, SpeciesContaminationCheck, GenotypeValidation,
                       SampleContaminationCheck, GenderValidation, CoverageStats)

    def _run(self):
        dir_with_linked_files = _link_results_files(self.dataset_name, self.job_dir, 'bcbio')
        return _output_data(self.dataset, self.job_dir, self.dataset_name, dir_with_linked_files)
