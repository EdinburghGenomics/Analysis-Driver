__author__ = 'mwham'
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
        self.crawler = report_generation.RunCrawler(
            self.run_id,
            SampleSheet(os.path.join(self.test_data, 'SampleSheet_analysis_driver.csv')),
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
        'properly_mapped_reads': 949154225,
        'duplicate_reads': 171911966,
        'sample_id': 'test_sample',
        'median_coverage': 30.156,
        'user_sample_id': 'test_sample',
        'pc_callable': 0.24392084973311048,
        'project_id': 'test_project',
        'bam_file_reads': 988805087,
        'mapped_reads': 975587288,
        'called_gender': 'male',
        'provided_gender': 'female',
        "species_contamination": {"contaminant_unique_mapped": {"Bos taurus": 1, "Felis catus": 4, "Gallus gallus": 1, "Mus musculus": 4, "Ovis aries": 2},
                                  "percent_unmapped": 1.06,
                                  "percent_unmapped_focal": 1.09,
                                  "total_reads_mapped": 100000
                                  },
        'coverage': {'mean': 438, 'median': 478, 'std_dev': 189}
    }

    def setUp(self):
        with patch('analysis_driver.report_generation.report_crawlers.get_user_sample_name',
                   return_value='test_sample'):
            with patch('analysis_driver.report_generation.report_crawlers.get_sex_from_lims',
                       return_value='female'):
                with patch('analysis_driver.reader.demultiplexing_parsers.get_species_from_sample',
                           return_value='Homo sapiens'):
                    self.crawler = report_generation.SampleCrawler('test_sample', 'test_project', self.test_data)

    def test_sample(self):
        self.compare_jsons(self.crawler.sample, self.expected_sample)