import os
import json
from analysis_driver.reader import demultiplexing_parsers as dm
from tests.test_analysisdriver import TestAnalysisDriver
from egcg_core.constants import ELEMENT_CONTAMINANT_UNIQUE_MAP, ELEMENT_PCNT_UNMAPPED_FOCAL,\
    ELEMENT_PCNT_UNMAPPED, ELEMENT_TOTAL_READS_MAPPED


class TestDemultiplexingStats(TestAnalysisDriver):
    def test(self):
        # testing multiplexed run stats parsing (previously from ConversionStats.xml and AdapterTrimming.txt')
        json_stat = os.path.join(self.assets_path, 'test_crawlers', 'Stats.json')
        with open(json_stat, 'r') as stats:
            json_data = json.load(stats)

        expected_barcode_per_lane = [
            ('LP6002014-DTP_A01', '1', 'ATTACTCG', 537099, 537099, 80564850, 72789430, 58579087),
            ('LP6002014-DTP_A02', '1', 'TCCGGAGA', 539999, 539999, 80999850, 73257234, 57184580),
            ('Undetermined', '1', 'unknown', 1696088, 80740, 12111000, 9937085, 7746149),
            ('LP6002014-DTP_A01', '2', 'ATTACTCG', 466412, 466412, 69961800, 61852875, 49303565),
            ('LP6002014-DTP_A02', '2', 'TCCGGAGA', 469731, 469731, 70459650, 62355216, 48408701),
            ('Undetermined', '2', 'unknown', 2335276, 122074, 18311100, 15157909, 12445247)
        ]
        expected_unknown_barcodes_per_lanes = [
            ('1', 'CCCCCCCC', '4820'),
            ('1', 'TGAAGCTA', '660'),
            ('2', 'CCCCCCCC', '21320'),
            ('2', 'NCCCCCCC', '9840')
        ]
        barcodes, unknowns, adapter_trimmed_by_id = dm.parse_json_stats(json_data, 'test_run_id')
        self.assertCountEqual(barcodes, expected_barcode_per_lane)
        self.assertCountEqual(unknowns, expected_unknown_barcodes_per_lanes)

        self.assertEqual(adapter_trimmed_by_id, {
            'test_run_id_1_ATTACTCG': {'read_1_trimmed_bases': 1056218, 'read_2_trimmed_bases': 1080753},
            'test_run_id_1_TCCGGAGA': {'read_1_trimmed_bases': 974374, 'read_2_trimmed_bases': 992558},
            'test_run_id_1_unknown': {'read_1_trimmed_bases': 100051, 'read_2_trimmed_bases': 122287},
            'test_run_id_2_ATTACTCG': {'read_1_trimmed_bases': 866438, 'read_2_trimmed_bases': 914963},
            'test_run_id_2_TCCGGAGA': {'read_1_trimmed_bases': 797469, 'read_2_trimmed_bases': 843785},
            'test_run_id_2_unknown': {'read_1_trimmed_bases': 112371, 'read_2_trimmed_bases': 143299}
        })

        # testing barcodeless run stats parsing (previously from ConversionStats.xml and AdapterTrimming.txt')
        json_stat = os.path.join(self.assets_path, 'test_crawlers', 'Stats_barcodeless.json')
        with open(json_stat, 'r') as stats:
            json_data = json.load(stats)

        expected_barcode_per_lane = [
            ('LP3645274-NTP_E01', '1', None, 621211104, 484863400, 72729510000, 65578821723, 46299068025),
            ('LP3645554-NTP_A08', '2', None, 621211104, 485960674, 72894101100, 66487334397, 49457604977)
        ]
        expected_unknown_barcodes_per_lanes = [
            ('1', 'unknown', '484863320'),
            ('2', 'unknown', '485960620')
        ]
        barcodes, unknowns, adapter_trimmed_by_id = dm.parse_json_stats(json_data, 'test_run_id')
        self.assertEqual(barcodes, expected_barcode_per_lane)
        self.assertEqual(unknowns, expected_unknown_barcodes_per_lanes)

        self.assertEqual(adapter_trimmed_by_id, {
            'test_run_id_1': {'read_1_trimmed_bases': 358146231, 'read_2_trimmed_bases': 389554370},
            'test_run_id_2': {'read_1_trimmed_bases': 735649209, 'read_2_trimmed_bases': 747065613}
        })

    def test_parse_seqtk_fqchk(self):
        fqchk_file = os.path.join(self.assets_path, '10015ATpool01_S1_L001_R1_001.fastq.gz.fqchk')
        nb_read, nb_base, lo_q, hi_q = dm.parse_seqtk_fqchk_file(fqchk_file, q_threshold=30)
        assert nb_read == 561151
        assert nb_base == 83750569
        assert lo_q == 8551190
        assert hi_q == 75199379

    def test_parse_fastqscreen_file1(self):
        screen_file = os.path.join(self.assets_path, 'test_sample_R1_screen.txt')
        result = dm.parse_fastqscreen_file(screen_file, 'Homo sapiens')
        assert result == {
            ELEMENT_PCNT_UNMAPPED: 1.06,
            ELEMENT_PCNT_UNMAPPED_FOCAL: 1.09,
            ELEMENT_TOTAL_READS_MAPPED: 100000,
            ELEMENT_CONTAMINANT_UNIQUE_MAP: {
                'Gallus gallus': 1, 'Felis catus': 4, 'Bos taurus': 1,
                'Ovis aries': 2, 'Mus musculus': 4, 'Homo sapiens': 74144
            }
        }

    def test_parse_fastqscreen_file2(self):
        screen_file = os.path.join(self.assets_path, 'test_sample_R1_screen.txt')
        result = dm.parse_fastqscreen_file(screen_file, 'Mellivora capensis')
        assert result == {
            'percent_unmapped': 1.06,
            'percent_unmapped_focal': 100.0,
            'total_reads_mapped': 100000,
            'contaminant_unique_mapped': {
                'Bos taurus': 1, 'Felis catus': 4, 'Gallus gallus': 1,
                'Homo sapiens': 74144, 'Ovis aries': 2, 'Mus musculus': 4
            }
        }

    def test_calculate_mean(self):
        hist_file = os.path.join(self.assets_path, 'test_sample.depth')
        hist = dm.collapse_histograms(dm.read_histogram_file(hist_file))
        test_mean = dm.calculate_mean(hist)
        assert test_mean == 438.8514851485148

    def test_calculate_median(self):
        hist_file = os.path.join(self.assets_path, 'test_sample.depth')
        hist = dm.collapse_histograms(dm.read_histogram_file(hist_file))
        test_median = dm.calculate_median(hist)
        assert test_median == 478

    def test_calculate_sd(self):
        hist_file = os.path.join(self.assets_path, 'test_sample.depth')
        hist = dm.collapse_histograms(dm.read_histogram_file(hist_file))
        test_sd = dm.calculate_sd(hist)
        assert test_sd == 189.1911391390011

    def test_calculate_bases_at_coverage(self):
        histogram = {5: 3, 10: 6, 15: 24, 20: 30, 25: 21, 30: 43, 35: 63}
        bases_5x, bases_15x, bases_30x = dm.calculate_bases_at_coverage(histogram)
        assert bases_5x == 187
        assert bases_15x == 157
        assert bases_30x == 63
        histogram = {5: 3, 10: 6, 15: 24, 20: 30, 25: 21}
        bases_5x, bases_15x, bases_30x = dm.calculate_bases_at_coverage(histogram)
        assert bases_30x == 0

    def test_get_n_percentile(self):
        histogram = {1: 5, 2: 2, 3: 4, 4: 6, 5: 3}
        assert dm.get_percentile(histogram, 50) == 3
        histogram = {1: 1, 2: 3, 3: 4, 4: 6, 5: 6}
        assert dm.get_percentile(histogram, 50) == 4
        histogram = {1: 3, 2: 3, 3: 4, 4: 6, 5: 4}
        assert dm.get_percentile(histogram, 50) == 3.5

    def test_get_coverage_statistics(self):
        hist_file = os.path.join(self.assets_path, 'test_sample.depth')
        mean, median, sd, coverage_pcs, cov, genome_size, evenness = dm.get_coverage_statistics(hist_file)
        assert mean == 438.8514851485148
        assert median == 478
        assert sd == 189.1911391390011
        assert genome_size == 101
        assert evenness == 0.8139335573648481
        assert coverage_pcs == {'percentile_5': 102, 'percentile_25': 279, 'percentile_50': 478,
                                'percentile_75': 625, 'percentile_95': 648}

    def test_calculate_size_genome(self):
        histogram = {1: 5, 2: 2, 3: 4, 4: 6, 5: 3}
        assert dm.calculate_size_genome(histogram) == 20

    def test_calculate_evenness(self):
        histogram = {1: 5, 2: 2, 3: 4, 4: 6, 5: 3}
        e = dm.calculate_evenness(histogram)
        assert e == 0.8
        histogram = {1: 5, 2: 2, 3: 4, 4: 6, 5: 3, 100: 5}
        assert dm.calculate_evenness(histogram) < e
        histogram = {0: 5000, 1: 5, 2: 2, 3: 4, 4: 6, 5: 3}
        assert dm.calculate_evenness(histogram) == 0

    def test_get_coverage_Y_chrom(self):
        hist_file = os.path.join(self.assets_path, 'test_sample_chrY.depth')
        mean = dm.get_coverage_y_chrom(hist_file)
        assert mean == 234.275

    def test_parse_welldup_file(self):
        welldup_file = os.path.join(self.assets_path, 'test_crawlers', 'test_run.well_dup')
        dup_per_lane = dm.parse_welldup_file(welldup_file)
        assert dup_per_lane == {1: 11.747, 2: 14.576, 4: 20.496, 5: 5.981, 6: 10.917, 7: 14.611, 8: 26.416}

    def test_parse_fastqFilterer_stats(self):
        fastqfilterer_stats = os.path.join(self.assets_path, 'RE_fastqfilterer.stats')
        expected_dict = {
            'r1i': 'RE_R1.fastq.gz',
            'r2i': 'RE_R2.fastq.gz',
            'r1o': 'RE_R1_filtered.fastq',
            'r2o': 'RE_R2_filtered.fastq',
            'read_pairs_checked': '20',
            'read_pairs_removed': '13',
            'read_pairs_remaining': '7',
            'remove_tiles': '1102,2202',
            'trim_r1': '14', 'trim_r2': '16'
        }
        assert dm.parse_fastq_filterer_stats(fastqfilterer_stats) == expected_dict

    def test_parse_interop_summary(self):
        interop_summary = os.path.join(self.assets_path, 'test_crawlers', 'test_run_dir', 'interop_summary.txt')

        res = dm.parse_interop_summary(interop_summary)
        expected_metrics = {
            'pc_clust_pf_r1': 74.36, 'pc_clust_pf_stdev_r1': 3.41, 'phasing_r1': 0.091, 'prephasing_r1': 0.06,
            'pc_q30_r1': 93.12, 'pc_aligned_r1': 0.51, 'pc_aligned_stdev_r1': 0.02, 'pc_error_r1': 0.38,
            'pc_error_stdev_r1': 0.13, 'intensity_c1_r1': 201.0, 'intensity_c1_stdev_r1': 7.0, 'pc_clust_pf_r2': 74.36,
            'pc_clust_pf_stdev_r2': 3.41, 'phasing_r2': 0.158, 'prephasing_r2': 0.07, 'pc_q30_r2': 80.67,
            'pc_aligned_r2': 0.51, 'pc_aligned_stdev_r2': 0.03, 'pc_error_r2': 0.61, 'pc_error_stdev_r2': 0.28,
            'intensity_c1_r2': 158.0, 'intensity_c1_stdev_r2': 10.0
        }
        assert res['1'] == expected_metrics

        interop_summary = os.path.join(self.assets_path, 'interop_summary_empty.txt')
        assert dm.parse_interop_summary(interop_summary) == {}
