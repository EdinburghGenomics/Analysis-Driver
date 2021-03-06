from os.path import join

from egcg_core.constants import ELEMENT_MEAN_INSERT_SIZE, ELEMENT_STD_DEV_INSERT_SIZE, ELEMENT_MEDIAN_INSERT_SIZE, \
    ELEMENT_MEDIAN_ABS_DEV_INSERT_SIZE

from tests.test_analysisdriver import TestAnalysisDriver
from analysis_driver.reader import mapping_stats_parsers as mp


class TestMappingStats(TestAnalysisDriver):
    expected_headers = '\t'.join(
        ['Sample', 'NO_CALL_NO_CALL', 'NO_CALL_HOM_REF', 'NO_CALL_HET', 'NO_CALL_HOM_VAR',
         'NO_CALL_UNAVAILABLE', 'NO_CALL_MIXED', 'HOM_REF_NO_CALL', 'HOM_REF_HOM_REF', 'HOM_REF_HET',
         'HOM_REF_HOM_VAR', 'HOM_REF_UNAVAILABLE', 'HOM_REF_MIXED', 'HET_NO_CALL', 'HET_HOM_REF', 'HET_HET',
         'HET_HOM_VAR', 'HET_UNAVAILABLE', 'HET_MIXED', 'HOM_VAR_NO_CALL', 'HOM_VAR_HOM_REF', 'HOM_VAR_HET',
         'HOM_VAR_HOM_VAR', 'HOM_VAR_UNAVAILABLE', 'HOM_VAR_MIXED', 'UNAVAILABLE_NO_CALL',
         'UNAVAILABLE_HOM_REF', 'UNAVAILABLE_HET', 'UNAVAILABLE_HOM_VAR', 'UNAVAILABLE_UNAVAILABLE',
         'UNAVAILABLE_MIXED', 'MIXED_NO_CALL', 'MIXED_HOM_REF', 'MIXED_HET', 'MIXED_HOM_VAR',
         'MIXED_UNAVAILABLE', 'MIXED_MIXED', 'Mismatching_Alleles']
    )
    expected_lines = [
        '\t'.join(
            ['ALL', '3', '0', '0', '0', '3859', '0', '0', '9', '0', '0', '34303', '0', '0', '0', '13', '0',
             '127', '0', '1', '0', '0', '6', '14', '0', '0', '0', '0', '0', '97', '0', '0', '0', '0', '0',
             '0', '0', '0']
        ),
        '\t'.join(
            ['T00001P001A01', '3', '0', '0', '0', '3859', '0', '0', '9', '0', '0', '34303', '0', '0', '0',
             '13', '0', '127', '0', '1', '0', '0', '6', '14', '0', '0', '0', '0', '0', '97', '0', '0', '0',
             '0', '0', '0', '0', '0']
        )
    ]

    def setUp(self):
        self.geno_val_file = join(self.assets_path, 'sample_data', 'T00001P001A01-validation.txt')
        self.vbi_selfSM_file = join(self.assets_path, 'sample_data', 'test_sample-chr22-vbi.selfSM')
        self.samtools_stats = join(self.assets_path, 'test_crawlers', 'samtools_stats.txt')
        self.vcf_stats_file = join(self.assets_path, 'test_crawlers', 'test_sample.vcf.stats')
        self.picard_markdup_file = join(self.assets_path, 'test_crawlers', 'test_picard_markdup.metrics')
        self.picard_markdup_file2 = join(self.assets_path, 'test_crawlers', 'test_picard2_markdup.metrics')
        self.picard_insertsize_file = join(self.assets_path, 'test_crawlers', 'test_picard_insertsize.metrics')

    def test_parse_genotype_concordance(self):
        table_type, headers, lines = mp.parse_genotype_concordance(self.geno_val_file)
        assert self.expected_headers == '\t'.join(headers.split())
        assert self.expected_lines == ['\t'.join(l.split()) for l in lines]

    def test_aggregate_genotype_concordance_lines(self):
        samples = mp.aggregate_genotype_concordance(self.expected_headers, self.expected_lines)
        assert set(samples.keys()) == {'ALL', 'T00001P001A01'}
        assert samples['T00001P001A01'] == {
            'matching_snps': 28, 'no_call_chip': 4, 'no_call_seq': 0, 'mismatching_snps': 0
        }

    def test_parse_vbi_selfSM(self):
        freemix = mp.parse_vbi_self_sm(self.vbi_selfSM_file)
        assert freemix == 0.00605

    def test_samtools_stats_parser(self):
        total_reads, mapped_reads, duplicate_reads, proper_pairs = mp.parse_samtools_stats(self.samtools_stats)
        assert total_reads == 7928618
        assert mapped_reads == 7892452
        assert duplicate_reads == 676698
        assert proper_pairs == 7741548

    def test_parse_vcf_stats(self):
        ti_tv, het_hom = mp.parse_vcf_stats(self.vcf_stats_file)
        assert ti_tv == 2.01
        assert het_hom == 1.6

    def test_picard_mark_dup(self):
        mapped_reads, duplicate_reads, opt_duplicate_reads, est_library_size = mp.parse_picard_mark_dup_metrics(self.picard_markdup_file)
        assert mapped_reads == 1059293
        assert duplicate_reads == 158975
        assert opt_duplicate_reads == 115966
        assert est_library_size == 5054793

        mapped_reads, duplicate_reads, opt_duplicate_reads, est_library_size = mp.parse_picard_mark_dup_metrics(self.picard_markdup_file2)
        assert mapped_reads == 1059293
        assert duplicate_reads == 158975
        assert opt_duplicate_reads == 115966
        assert est_library_size is None

    def test_parse_picard_insert_size_metrics(self):
        insert_types = mp.parse_picard_insert_size_metrics(self.picard_insertsize_file)
        assert 'FR' in insert_types
        assert insert_types['FR'][ELEMENT_MEAN_INSERT_SIZE] == 446.299806
        assert insert_types['FR'][ELEMENT_STD_DEV_INSERT_SIZE] == 113.966554
        assert insert_types['FR'][ELEMENT_MEDIAN_INSERT_SIZE] == 439
        assert insert_types['FR'][ELEMENT_MEDIAN_ABS_DEV_INSERT_SIZE] == 69
