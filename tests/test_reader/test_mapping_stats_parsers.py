import os
from analysis_driver.reader.mapping_stats_parsers import parse_genotype_concordance
from tests.test_analysisdriver import TestAnalysisDriver

__author__ = 'tcezard'


class TestMappingStats(TestAnalysisDriver):
    def setUp(self):
        self.geno_val_file = os.path.join(self.assets_path, 'sample_data', 'T00001P001A01-validation.txt')

    def test_parse_genotype_concordance(self):
        expected_headers = '\t'.join(
                ["Sample", "NO_CALL_NO_CALL", "NO_CALL_HOM_REF", "NO_CALL_HET", "NO_CALL_HOM_VAR",
                   "NO_CALL_UNAVAILABLE", "NO_CALL_MIXED", "HOM_REF_NO_CALL", "HOM_REF_HOM_REF",
                   "HOM_REF_HET", "HOM_REF_HOM_VAR", "HOM_REF_UNAVAILABLE", "HOM_REF_MIXED", "HET_NO_CALL",
                   "HET_HOM_REF", "HET_HET", "HET_HOM_VAR", "HET_UNAVAILABLE", "HET_MIXED", "HOM_VAR_NO_CALL",
                   "HOM_VAR_HOM_REF", "HOM_VAR_HET", "HOM_VAR_HOM_VAR", "HOM_VAR_UNAVAILABLE", "HOM_VAR_MIXED",
                   "UNAVAILABLE_NO_CALL", "UNAVAILABLE_HOM_REF", "UNAVAILABLE_HET", "UNAVAILABLE_HOM_VAR",
                   "UNAVAILABLE_UNAVAILABLE", "UNAVAILABLE_MIXED", "MIXED_NO_CALL", "MIXED_HOM_REF", "MIXED_HET",
                   "MIXED_HOM_VAR", "MIXED_UNAVAILABLE", "MIXED_MIXED", "Mismatching_Alleles"]
        )
        value1 = '\t'.join(
                ["ALL", "3", "0", "0", "0", "3859", "0", "0", "9", "0", "0", "34303", "0", "0", "0", "13",
                 "0", "127", "0", "1", "0", "0", "6", "14", "0", "0", "0", "0", "0", "97", "0", "0", "0",
                 "0", "0", "0", "0", "0"]
                           )
        value2 = '\t'.join(["T00001P001A01", "3", "0", "0", "0", "3859", "0", "0", "9", "0", "0", "34303", "0", "0",
                            "0", "13", "0", "127", "0", "1", "0", "0", "6", "14", "0", "0", "0", "0", "0", "97", "0",
                            "0", "0", "0", "0", "0", "0", "0"])
        expected_lines = [value1, value2]

        table_type, headers, lines = parse_genotype_concordance(self.geno_val_file)
        assert expected_headers == '\t'.join(headers.split())
        assert expected_lines == ['\t'.join(l.split()) for l in lines]


    def test_aggregate_genotype_concordance_lines(self):
        headers = '\t'.join(
                ["Sample NO_CALL_NO_CALL", "NO_CALL_HOM_REF", "NO_CALL_HET", "NO_CALL_HOM_VAR",
                   "NO_CALL_UNAVAILABLE", "NO_CALL_MIXED", "HOM_REF_NO_CALL", "HOM_REF_HOM_REF",
                   "HOM_REF_HET", "HOM_REF_HOM_VAR", "HOM_REF_UNAVAILABLE", "HOM_REF_MIXED", "HET_NO_CALL",
                   "HET_HOM_REF", "HET_HET", "HET_HOM_VAR", "HET_UNAVAILABLE", "HET_MIXED", "HOM_VAR_NO_CALL",
                   "HOM_VAR_HOM_REF", "HOM_VAR_HET", "HOM_VAR_HOM_VAR", "HOM_VAR_UNAVAILABLE", "HOM_VAR_MIXED",
                   "UNAVAILABLE_NO_CALL", "UNAVAILABLE_HOM_REF", "UNAVAILABLE_HET", "UNAVAILABLE_HOM_VAR",
                   "UNAVAILABLE_UNAVAILABLE", "UNAVAILABLE_MIXED", "MIXED_NO_CALL", "MIXED_HOM_REF", "MIXED_HET",
                   "MIXED_HOM_VAR", "MIXED_UNAVAILABLE", "MIXED_MIXED", "Mismatching_Alleles"]
        )
        value1 = '\t'.join(
                ["ALL", "3", "0", "0", "0", "3859", "0", "0", "9", "0", "0", "34303", "0", "0", "0", "13",
                 "0", "127", "0", "1", "0", "0", "6", "14", "0", "0", "0", "0", "0", "97", "0", "0", "0",
                 "0", "0", "0", "0", "0"]
                           )
        value2 = '\t'.join(["T00001P001A01", "3", "0", "0", "0", "3859", "0", "0", "9", "0", "0", "34303", "0", "0",
                            "0", "13", "0", "127", "0", "1", "0", "0", "6", "14", "0", "0", "0", "0", "0", "97", "0",
                            "0", "0", "0", "0", "0", "0", "0"])
        lines = [value1, value2]
        samples = aggregate_genotype_concordance(headers, lines)
        assert set(samples.keys()) == set(['ALL', 'T00001P001A01'])
        assert dict(samples['T00001P001A01']) == {
            'matching_snps': 28, 'no_call_chip': 4, 'no_call_seq': 0, 'mismatching_snps': 0
        }
