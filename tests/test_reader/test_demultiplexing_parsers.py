import os
from analysis_driver.reader.demultiplexing_parsers import parse_seqtk_fqchk_file
from tests.test_analysisdriver import TestAnalysisDriver

__author__ = 'tcezard'


class Test_demultiplexing_stats(TestAnalysisDriver):

    def test_parse_genotype_concordance(self):
        fqchk_file = os.path.join(self.assets_path, '10015ATpool01_S1_L001_R1_001.fastq.gz.fqchk')
        results = parse_seqtk_fqchk_file(fqchk_file, q_threshold=30)
        print(results)



