import os
from analysis_driver.reader.demultiplexing_parsers import parse_seqtk_fqchk_file
from analysis_driver.reader.demultiplexing_parsers import parse_fastqscreen_file
from analysis_driver.reader.demultiplexing_parsers import get_fastqscreen_results
from tests.test_analysisdriver import TestAnalysisDriver
from analysis_driver.constants import ELEMENT_CONTAMINANT_UNIQUE_MAP, ELEMENT_PCNT_UNMAPPED_FOCAL, ELEMENT_PCNT_UNMAPPED
from unittest.mock import patch
__author__ = 'tcezard'


class TestDemultiplexingStats(TestAnalysisDriver):

    def test_parse_seqtk_fqchk(self):
        fqchk_file = os.path.join(self.assets_path, '10015ATpool01_S1_L001_R1_001.fastq.gz.fqchk')
        nb_read, nb_base, lo_q, hi_q = parse_seqtk_fqchk_file(fqchk_file, q_threshold=30)
        assert nb_read == 561151
        assert nb_base == 83750569
        assert lo_q == 8551190
        assert hi_q == 75199379

    def test_parse_fastqscreen_file1(self):
        testFile = os.path.join(self.assets_path, "fastqscreenTestOutput.txt")
        result = parse_fastqscreen_file(testFile, 'Homo sapiens')
        assert result == {ELEMENT_PCNT_UNMAPPED: 1.06, ELEMENT_PCNT_UNMAPPED_FOCAL: 1.09, ELEMENT_CONTAMINANT_UNIQUE_MAP: {'Gallus gallus': 1, 'Felis catus': 4, 'Bos taurus': 1, 'Ovis aries': 2, 'Mus musculus': 4}}

    def test_parse_fastqscreen_file2(self):
        testFile = os.path.join(self.assets_path, "fastqscreenTestOutput.txt")
        result = parse_fastqscreen_file(testFile, 'Mellivora capensis')
        assert result == [100, 100, 100]

    @patch('analysis_driver.reader.demultiplexing_parsers.get_species_from_sample', autospec=True)
    def test_get_fastqscreen_results1(self, mocked_species_sample):
        testFile = os.path.join(self.assets_path, "fastqscreenTestOutput.txt")
        mocked_species_sample.return_value = None
        result = get_fastqscreen_results(testFile, 'testSampleID')
        assert result == [100, 100, 100]

    @patch('analysis_driver.reader.demultiplexing_parsers.get_species_from_sample', autospec=True)
    def test_get_fastqscreen_results2(self, mocked_species_sample):
        testFile = os.path.join(self.assets_path, "fastqscreenTestOutput.txt")
        mocked_species_sample.return_value = 'Homo sapiens'
        result = get_fastqscreen_results(testFile, 'testSampleID')
        assert result == {ELEMENT_PCNT_UNMAPPED: 1.06, ELEMENT_PCNT_UNMAPPED_FOCAL: 1.09, ELEMENT_CONTAMINANT_UNIQUE_MAP: {'Gallus gallus': 1, 'Felis catus': 4, 'Bos taurus': 1, 'Ovis aries': 2, 'Mus musculus': 4}}

