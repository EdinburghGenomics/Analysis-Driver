import os.path
import json
from unittest.mock import patch
from tests.test_analysisdriver import TestAnalysisDriver
from analysis_driver import report_generation
from analysis_driver.reader import SampleSheet


class TestCrawler(TestAnalysisDriver):
    @property
    def test_data(self):
        return os.path.join(self.assets_path, 'test_crawlers')

    @classmethod
    def compare_jsons(cls, observed, expected):
        cls._sort_lists(observed)
        cls._sort_lists(expected)
        o = json.dumps(observed, sort_keys=True)
        e = json.dumps(expected, sort_keys=True)
        if o != e:
            print('observed:')
            print(o)
            print('expected:')
            print(e)
            raise AssertionError

    @classmethod
    def _sort_lists(cls, d):
        for k, v in d.items():
            if type(v) is list:
                d[k] = sorted(v)
            elif type(v) is dict:
                cls._sort_lists(v)

class TestRunCrawler(TestCrawler):
    run_id = '150723_E00306_0025_BHCHK3CCXX'
    _expected_output = None

    @property
    def expected_output(self):
        if self._expected_output is None:
            self._expected_output = json.load(
                open(os.path.join(self.test_data, 'expected_run_crawler_data.json'))
            )
        return self._expected_output

    def setUp(self):
        with patch('analysis_driver.reader.demultiplexing_parsers.run_sample_lane_to_barcode',
                   return_value={"150723_E00306_0025_BHCHK3CCXX_1_unknown": {'read_1_trimmed_bases': 184380158, 'read_2_trimmed_bases': 172552099},
                                 "150723_E00306_0025_BHCHK3CCXX_2_unknown": {'read_1_trimmed_bases': 48149799, 'read_2_trimmed_bases': 48818739},
                                 "150723_E00306_0025_BHCHK3CCXX_1_TCCGGAGA": {'read_1_trimmed_bases': 1088149481, 'read_2_trimmed_bases': 1034179505},
                                 "150723_E00306_0025_BHCHK3CCXX_2_TCCGGAGA": {'read_1_trimmed_bases': 398993728, 'read_2_trimmed_bases': 391621660},
                                 "150723_E00306_0025_BHCHK3CCXX_1_ATTACTCG": {'read_1_trimmed_bases': 714309214, 'read_2_trimmed_bases': 684692293},
                                 "150723_E00306_0025_BHCHK3CCXX_2_ATTACTCG": {'read_1_trimmed_bases': 284712861, 'read_2_trimmed_bases': 282625840}}):

            self.crawler = report_generation.RunCrawler(
                self.run_id,
                SampleSheet(os.path.join(self.test_data, 'SampleSheet_analysis_driver.csv')),
                os.path.join(self.test_data, 'AdapterTrimming.txt'),
                os.path.join(self.test_data, 'ConversionStats.xml')
            )

    def test_barcodes_info(self):
        self.compare_jsons(dict(self.crawler.barcodes_info), self.expected_output['barcodes_info'])

    def test_unexpected_barcodes(self):
        self.compare_jsons(dict(self.crawler.unexpected_barcodes), self.expected_output['unexpected_barcodes'])

    def test_libraries(self):
        self.compare_jsons(dict(self.crawler.libraries), self.expected_output['libraries'])

    def test_lanes(self):
        self.compare_jsons(dict(self.crawler.lanes), self.expected_output['lanes'])

    def test_run(self):
        self.compare_jsons(dict(self.crawler.run), self.expected_output['run'])

    def test_projects(self):
        self.compare_jsons(dict(self.crawler.projects), self.expected_output['projects'])

class TestSampleCrawler(TestCrawler):
    expected_sample = {
        'properly_mapped_reads': 7741548,
        'duplicate_reads': 676698,
        'sample_id': 'test_sample',
        'median_coverage': 478,
        'user_sample_id': 'test_sample',
        'pc_callable': 0.24392084973311048,
        'project_id': 'test_project',
        'bam_file_reads': 7928618,
        'mapped_reads': 7892452,
        'called_gender': 'male',
        'provided_gender': 'female',
        'species_contamination': {
            'contaminant_unique_mapped': {
                'Bos taurus': 1,
                'Felis catus': 4,
                'Gallus gallus': 1,
                'Mus musculus': 4,
                'Ovis aries': 2,
                'Homo sapiens': 74144
            },
            'percent_unmapped': 1.06,
            'percent_unmapped_focal': 1.09,
            'total_reads_mapped': 100000
        },
        "sample_contamination": {"het_hom_ratio": 1.6, "ti_tv_ratio": 2.01},
        'gender_validation': {'hetX': '0.10'},
        'coverage': {'median': 478,
                     'std_dev': 189.1911391390011,
                     'coverage_percentiles': {"percentile_25": 279, "percentile_5": 102, "percentile_50": 478, "percentile_75": 625, "percentile_95": 648},
                     'mean': 438.8514851485148,
                     'bases_at_coverage': {'bases_at_5X': 100, 'bases_at_30X': 99, 'bases_at_15X': 100}}
    }

    def setUp(self):
        with patch('analysis_driver.report_generation.report_crawlers.get_user_sample_name',
                   return_value='test_sample'):
            with patch('analysis_driver.report_generation.report_crawlers.get_sample_gender',
                       return_value='female'):
                with patch('analysis_driver.reader.demultiplexing_parsers.get_species_from_sample',
                           return_value='Homo sapiens'):
                    self.crawler = report_generation.SampleCrawler('test_sample', 'test_project', self.test_data)

    def test_sample(self):
        self.compare_jsons(self.crawler.sample, self.expected_sample)
