import os
from analysis_driver.reader.demultiplexing_parsers import parse_seqtk_fqchk_file
from tests.test_analysisdriver import TestAnalysisDriver

__author__ = 'tcezard'


class TestDemultiplexingStats(TestAnalysisDriver):

    def test_parse_seqtk_fqchk(self):
        fqchk_file = os.path.join(self.assets_path, '10015ATpool01_S1_L001_R1_001.fastq.gz.fqchk')
        nb_read, nb_base, lo_q, hi_q = parse_seqtk_fqchk_file(fqchk_file, q_threshold=30)
        assert nb_read == 561151
        assert nb_base == 83750569
        assert lo_q == 8551190
        assert hi_q == 75199379
