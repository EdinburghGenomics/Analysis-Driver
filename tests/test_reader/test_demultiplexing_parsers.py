import os
from analysis_driver.reader.demultiplexing_parsers import parse_seqtk_fqchk_file, parse_conversion_stats
from analysis_driver.reader.demultiplexing_parsers import parse_fastqscreen_file
from analysis_driver.reader.demultiplexing_parsers import get_fastqscreen_results
from analysis_driver.reader.demultiplexing_parsers import calculate_mean, calculate_median, calculate_sd, get_coverage_statistics
from tests.test_analysisdriver import TestAnalysisDriver
from analysis_driver.constants import ELEMENT_CONTAMINANT_UNIQUE_MAP, ELEMENT_PCNT_UNMAPPED_FOCAL, ELEMENT_PCNT_UNMAPPED, ELEMENT_TOTAL_READS_MAPPED
from unittest.mock import patch
__author__ = 'tcezard'


class TestDemultiplexingStats(TestAnalysisDriver):

    def test(self):
        conversion_stat = os.path.join(self.assets_path, 'test_crawlers', 'ConversionStats.xml')
        expected_barcode_per_lane = [
            ('default', 'Undetermined', '1', 'unknown', 1696088, 80740, 12111000, 9937085, 7746149),
            ('default', 'Undetermined', '2', 'unknown', 2335276, 122074, 18311100, 15157909, 12445247),
            ('10015AT', 'LP6002014-DTP_A01', '1', 'ATTACTCG', 537099, 537099, 80564850, 72789430, 58579087),
            ('10015AT', 'LP6002014-DTP_A01', '2', 'ATTACTCG', 466412, 466412, 69961800, 61852875, 49303565),
            ('10015AT', 'LP6002014-DTP_A02', '1', 'TCCGGAGA', 539999, 539999, 80999850, 73257234, 57184580),
            ('10015AT', 'LP6002014-DTP_A02', '2', 'TCCGGAGA', 469731, 469731, 70459650, 62355216, 48408701)
        ]
        expected_unknown_barcodes_per_lanes = [
            ('1', 'CCCCCCCC', '4695'),
            ('1', 'TGAAGCTA', '828'),
            ('2', 'CCCCCCCC', '22153'),
            ('2', 'NCCCCCCC', '9133')
        ]
        barcodes_per_lane, top_unknown_barcodes_per_lanes, barcodeless_per_lane = parse_conversion_stats(conversion_stat, has_barcode=True)
        assert barcodes_per_lane == expected_barcode_per_lane
        assert top_unknown_barcodes_per_lanes == expected_unknown_barcodes_per_lanes
        assert barcodeless_per_lane == []

        conversion_stat = os.path.join(self.assets_path, 'test_crawlers', 'ConversionStats_barcodeless.xml')
        expected_barcodeless_per_lane = [
            ('Corriell_2016-02-03', 'LP0000038-NTP_C01', '1', 'all', 6470949, 3794190, 569128500, 481721486, 344363516),
            ('Corriell_2016-02-03', 'LP0000038-NTP_C02', '2', 'all', 6470949, 4254880, 638232000, 542343434, 364820834)
        ]
        barcodes_per_lane, top_unknown_barcodes_per_lanes, barcodeless_per_lane = parse_conversion_stats(conversion_stat, has_barcode=False)
        assert barcodeless_per_lane == expected_barcodeless_per_lane



    def test_parse_seqtk_fqchk(self):
        fqchk_file = os.path.join(self.assets_path, '10015ATpool01_S1_L001_R1_001.fastq.gz.fqchk')
        nb_read, nb_base, lo_q, hi_q = parse_seqtk_fqchk_file(fqchk_file, q_threshold=30)
        assert nb_read == 561151
        assert nb_base == 83750569
        assert lo_q == 8551190
        assert hi_q == 75199379

    def test_parse_fastqscreen_file1(self):
        testFile = os.path.join(self.assets_path, "test_sample_R1_screen.txt")
        result = parse_fastqscreen_file(testFile, 'Homo sapiens')
        assert result == {ELEMENT_PCNT_UNMAPPED: 1.06, ELEMENT_PCNT_UNMAPPED_FOCAL: 1.09, ELEMENT_TOTAL_READS_MAPPED: 100000, ELEMENT_CONTAMINANT_UNIQUE_MAP: {'Gallus gallus': 1, 'Felis catus': 4, 'Bos taurus': 1, 'Ovis aries': 2, 'Mus musculus': 4, 'Homo sapiens': 74144}}

    def test_parse_fastqscreen_file2(self):
        testFile = os.path.join(self.assets_path, "test_sample_R1_screen.txt")
        result = parse_fastqscreen_file(testFile, 'Mellivora capensis')
        assert result is None

    @patch('analysis_driver.reader.demultiplexing_parsers.get_species_from_sample', autospec=True)
    def test_get_fastqscreen_results1(self, mocked_species_sample):
        testFile = os.path.join(self.assets_path, "test_sample_R1_screen.txt")
        mocked_species_sample.return_value = None
        result = get_fastqscreen_results(testFile, 'testSampleID')
        assert result is None

    @patch('analysis_driver.reader.demultiplexing_parsers.get_species_from_sample', autospec=True)
    def test_get_fastqscreen_results2(self, mocked_species_sample):
        testFile = os.path.join(self.assets_path, "test_sample_R1_screen.txt")
        mocked_species_sample.return_value = 'Homo sapiens'
        result = get_fastqscreen_results(testFile, 'testSampleID')
        assert result == {ELEMENT_PCNT_UNMAPPED: 1.06, ELEMENT_TOTAL_READS_MAPPED: 100000, ELEMENT_PCNT_UNMAPPED_FOCAL: 1.09, ELEMENT_CONTAMINANT_UNIQUE_MAP: {'Gallus gallus': 1, 'Felis catus': 4, 'Bos taurus': 1, 'Ovis aries': 2, 'Mus musculus': 4, 'Homo sapiens': 74144}}

    def test_calculate_mean(self):
        hist_file = os.path.join(self.assets_path, 'test_sample.hist')
        test_mean = calculate_mean(hist_file)
        assert test_mean == 438

    def test_calculate_median(self):
        hist_file = os.path.join(self.assets_path, 'test_sample.hist')
        test_median = calculate_median(hist_file)
        assert test_median == 478

    def test_calculate_sd(self):
        hist_file = os.path.join(self.assets_path, 'test_sample.hist')
        test_sd = calculate_sd(hist_file)
        assert test_sd == 189

    def test_get_coverage_statistics(self):
        hist_file = os.path.join(self.assets_path, 'test_sample.hist')
        mean, median, sd = get_coverage_statistics(hist_file)
        assert mean == 438
        assert median == 478
        assert sd == 189
