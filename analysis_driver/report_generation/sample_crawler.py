from egcg_core import util, clarity
from egcg_core.constants import *  # pylint: disable=unused-import
from egcg_core.rest_communication import post_or_patch as pp
from analysis_driver.reader import demultiplexing_parsers as dm, mapping_stats_parsers as mp
from analysis_driver.config import output_file_config
from .crawler import Crawler, sex_alias


class SampleCrawler(Crawler):
    def __init__(self, sample_id,  project_id, data_dir, output_fileset, post_pipeline=False):
        self.sample_id = sample_id
        self.user_sample_id = clarity.get_user_sample_name(sample_id, lenient=True)
        self.project_id = project_id
        self.all_info = []
        self.data_dir = data_dir
        self.post_pipeline = post_pipeline
        self.output_fileset = output_fileset
        self.sample = self._populate_lib_info()

    def get_output_file(self, outfile_id):
        if self.post_pipeline:
            fp = output_file_config.output_dir_file(self.output_fileset, outfile_id)
        else:
            fp = output_file_config.job_dir_file(self.output_fileset, outfile_id)

        self.debug('Searching for %s in %s/%s' % (outfile_id, self.data_dir, fp))
        if fp:
            return util.find_file(
                self.data_dir,
                fp.format(sample_id=self.sample_id, user_sample_id=self.user_sample_id)
            )

    def _populate_lib_info(self):
        sample = {
            ELEMENT_SAMPLE_INTERNAL_ID: self.sample_id,
            ELEMENT_PROJECT_ID: self.project_id,
        }
        sample.update(self.get_sample_information_from_lims(self.sample_id))

        samtools_path = self.get_output_file('samtools_stats')
        if samtools_path:
            (total_reads, mapped_reads,
             duplicate_reads, proper_pairs) = mp.parse_samtools_stats(samtools_path)

            sample[ELEMENT_NB_READS_IN_BAM] = total_reads
            sample[ELEMENT_NB_MAPPED_READS] = mapped_reads
            sample[ELEMENT_NB_DUPLICATE_READS] = duplicate_reads
            sample[ELEMENT_NB_PROPERLY_MAPPED] = proper_pairs
        else:
            self.critical('Missing samtools_stats.txt')

        bed_file_path = self.get_output_file('sort_callable')
        if bed_file_path:
            coverage_per_type = mp.parse_callable_bed_file(bed_file_path)
            callable_bases = coverage_per_type.get('CALLABLE')
            total = sum(coverage_per_type.values())
            sample[ELEMENT_PC_BASES_CALLABLE] = callable_bases/total
        else:
            self.critical('Missing *-sort-callable.bed')

        sex_file_path = self.get_output_file('sex_check')
        if sex_file_path:
            with open(sex_file_path) as f:
                sex, het_x = f.read().strip().split()
                sample[ELEMENT_SEX_VALIDATION][ELEMENT_CALLED_SEX] = sex_alias(sex)
                sample[ELEMENT_SEX_VALIDATION][ELEMENT_SEX_HETX] = het_x

        genotype_validation_path = self.get_output_file('genoval')
        if genotype_validation_path:
            genotyping_results = mp.parse_and_aggregate_genotype_concordance(genotype_validation_path)
            genotyping_result = genotyping_results.get(self.sample_id)
            if genotyping_result:
                sample[ELEMENT_GENOTYPE_VALIDATION] = genotyping_result
            else:
                self.critical('Sample %s not found in file %s', self.sample_id, genotype_validation_path)

        species_contamination_path = self.get_output_file('r1_fastqscreen')
        if species_contamination_path:
            species_contamination_result = dm.parse_fastqscreen_file(species_contamination_path,
                                                                     sample[ELEMENT_SAMPLE_SPECIES])
            if species_contamination_result:
                sample[ELEMENT_SPECIES_CONTAMINATION] = species_contamination_result
            else:
                self.critical('Contamination check unavailable for %s', self.sample_id)

        sample_contamination_path = self.get_output_file('self_sm')
        if sample_contamination_path:
            freemix = mp.parse_vbi_self_sm(sample_contamination_path)
            if freemix is not None:
                sample[ELEMENT_SAMPLE_CONTAMINATION] = {ELEMENT_FREEMIX: freemix}
            else:
                self.critical('freemix results from validateBamId are not available for %s', self.sample_id)

        coverage_statistics_path = self.get_output_file('depth_file')
        if coverage_statistics_path:
            mean, median, sd, coverage_percentiles, bases_at_coverage, \
             genome_size, evenness = dm.get_coverage_statistics(coverage_statistics_path)
            coverage_statistics = {
                ELEMENT_MEAN_COVERAGE: mean,
                ELEMENT_MEDIAN_COVERAGE_SAMTOOLS: median,
                ELEMENT_COVERAGE_SD: sd,
                ELEMENT_COVERAGE_PERCENTILES: coverage_percentiles,
                ELEMENT_BASES_AT_COVERAGE: bases_at_coverage,
                ELEMENT_SAMPLE_GENOME_SIZE: genome_size,
                ELEMENT_COVERAGE_EVENNESS: evenness
            }
            sample[ELEMENT_COVERAGE_STATISTICS] = coverage_statistics
            sample[ELEMENT_MEDIAN_COVERAGE] = median
            if ELEMENT_SEX_VALIDATION in sample:
                cov_y = dm.get_coverage_y_chrom(coverage_statistics_path)
                if cov_y:
                    sample[ELEMENT_SEX_VALIDATION][ELEMENT_SEX_COVY] = cov_y
        else:
            self.critical('coverage statistics unavailable for %s', self.sample_id)

        vcf_stats_path = self.get_output_file('vcf_stats')
        if vcf_stats_path:
            ti_tv, het_hom = mp.parse_vcf_stats(vcf_stats_path)
            if ELEMENT_SAMPLE_CONTAMINATION in sample:
                sample[ELEMENT_SAMPLE_CONTAMINATION][ELEMENT_SNPS_TI_TV] = ti_tv
                sample[ELEMENT_SAMPLE_CONTAMINATION][ELEMENT_SNPS_HET_HOM] = het_hom
            else:
                sample[ELEMENT_SAMPLE_CONTAMINATION] = {ELEMENT_SNPS_TI_TV: ti_tv, ELEMENT_SNPS_HET_HOM: het_hom}
        return sample

    def send_data(self):
        return pp('samples', [self.sample], ELEMENT_SAMPLE_INTERNAL_ID)
