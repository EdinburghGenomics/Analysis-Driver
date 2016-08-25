from os.path import join
from analysis_driver import quality_control as qc
from analysis_driver.pipeline import Stage
from analysis_driver.pipeline.common import MergeFastqs, Fastqc, CoverageStats
from analysis_driver.driver import _run_bcbio, _link_results_files, _output_data
from egcg_core import util


class BCBioStage(Stage):
    @property
    def fastq_pair(self):
        return util.find_files(self.job_dir, 'merged', self.dataset.user_sample_id + '_R?.fastq.gz')

    def _run(self):
        raise NotImplementedError


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


class DataOutput(BCBioStage):
    @property
    def previous_stages(self):
        return (
            Fastqc(previous_stages=MergeFastqs, fastqs=self.fastq_pair),
            SpeciesContaminationCheck,
            GenotypeValidation,
            SampleContaminationCheck,
            GenderValidation,
            CoverageStats(
                previous_stages=BCBio,
                bam_file=join(
                    'samples_%s-merged' % self.dataset_name,
                    'final',
                    self.dataset.user_sample_id,
                    self.dataset.user_sample_id + '-ready.bam'
                )
            )
        )

    def _run(self):
        dir_with_linked_files = _link_results_files(self.dataset_name, self.job_dir, 'bcbio')
        return _output_data(self.dataset, self.job_dir, self.dataset_name, dir_with_linked_files)
