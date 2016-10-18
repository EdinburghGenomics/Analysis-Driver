import os
from analysis_driver.reader.demultiplexing_parsers import parse_seqtk_fqchk_file, parse_conversion_stats, \
    parse_welldup_file, get_percentiles, read_histogram_file, collapse_histograms, get_coverage_Y_chrom
from analysis_driver.reader.demultiplexing_parsers import parse_fastqscreen_file
from analysis_driver.reader.demultiplexing_parsers import calculate_mean, calculate_median, calculate_sd, get_coverage_statistics, calculate_bases_at_coverage
from tests.test_analysisdriver import TestAnalysisDriver
from egcg_core.constants import ELEMENT_CONTAMINANT_UNIQUE_MAP, ELEMENT_PCNT_UNMAPPED_FOCAL, ELEMENT_PCNT_UNMAPPED, ELEMENT_TOTAL_READS_MAPPED
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

    def test_calculate_mean(self):
        hist_file = os.path.join(self.assets_path, 'test_sample.depth')
        hist = collapse_histograms(read_histogram_file(hist_file))
        test_mean = calculate_mean(hist)
        assert test_mean == 438.8514851485148

    def test_calculate_median(self):
        hist_file = os.path.join(self.assets_path, 'test_sample.depth')
        hist = collapse_histograms(read_histogram_file(hist_file))
        test_median = calculate_median(hist)
        assert test_median == 478

    def test_calculate_sd(self):
        hist_file = os.path.join(self.assets_path, 'test_sample.depth')
        hist = collapse_histograms(read_histogram_file(hist_file))
        test_sd = calculate_sd(hist)
        assert test_sd == 189.1911391390011


    def test_calculate_bases_at_coverage(self):
        histogram = {5:3, 10:6, 15:24, 20:30, 25:21, 30:43, 35:63}
        bases_5X, bases_15X, bases_30X = calculate_bases_at_coverage(histogram)
        assert bases_5X == 187
        assert bases_15X == 157
        assert bases_30X == 63
        histogram = {5:3, 10:6, 15:24, 20:30, 25:21}
        bases_5X, bases_15X, bases_30X = calculate_bases_at_coverage(histogram)
        assert bases_30X == 0

    def test_get_percentiles(self):
        histogram={1:5, 2:2, 3:4, 4:6, 5:3}
        assert get_percentiles(histogram, 50) == 3
        histogram={1:1, 2:3, 3:4, 4:6, 5:6}
        assert get_percentiles(histogram, 50) == 4
        histogram={1:3, 2:3, 3:4, 4:6, 5:4}
        assert get_percentiles(histogram, 50) == 3.5


    def test_get_coverage_statistics(self):
        hist_file = os.path.join(self.assets_path, 'test_sample.depth')
        mean, median, sd, coverage_percentiles, bases_at_coverage = get_coverage_statistics(hist_file)
        assert mean == 438.8514851485148
        assert median == 478
        assert sd == 189.1911391390011
        assert coverage_percentiles == {'percentile_5': 102, 'percentile_25': 279, 'percentile_50': 478, 'percentile_75': 625, 'percentile_95': 648}

    def test_get_coverage_Y_chrom(self):
        hist_file = os.path.join(self.assets_path, 'test_sample_chrY.depth')
        mean = get_coverage_Y_chrom(hist_file)
        assert mean == 234.275

    def test_parse_welldup_file(self):
        welldup_file = os.path.join(self.assets_path, 'test_crawlers', 'test_run.well_dup')
        dup_per_lane = parse_welldup_file(welldup_file)
        assert dup_per_lane == {1: 11.747, 2: 14.576, 3: 0, 4: 20.496, 5: 5.981, 6: 10.917, 7: 14.611, 8: 26.416}