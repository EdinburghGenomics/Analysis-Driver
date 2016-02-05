import os
from analysis_driver.reader.mapping_stats_parsers import parse_genotype_concordance
from tests.test_analysisdriver import TestAnalysisDriver

__author__ = 'tcezard'


class Test_mapping_stats(TestAnalysisDriver):
    def setUp(self):
        self.geno_val_file = os.path.join(self.assets_path,'sample_data','T00001P001A01-validation.txt')

    def test_parse_genotype_concordance(self):
        samples = parse_genotype_concordance(self.geno_val_file)
        assert set(samples.keys()) == set(['ALL', 'T00001P001A01'])
        assert dict(samples['T00001P001A01']) == {'matching_snps': 28, 'no_call_chip': 4, 'no_call_seq': 0, 'mismatching_snps': 0}




