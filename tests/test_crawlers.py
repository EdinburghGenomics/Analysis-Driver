import os.path
import json
from unittest.mock import patch
from egcg_core import constants as c
from tests.test_analysisdriver import TestAnalysisDriver, NamedMock
from analysis_driver import report_generation
from analysis_driver.config import OutputFileConfiguration

ppath = 'analysis_driver.report_generation.'


class TestCrawler(TestAnalysisDriver):
    @property
    def test_data(self):
        return os.path.join(self.assets_path, 'test_crawlers')

    @classmethod
    def compare_jsons(cls, observed, expected):
        cls._sort_lists(observed)
        cls._sort_lists(expected)
        o = json.dumps(observed, indent=4, sort_keys=True)
        e = json.dumps(expected, indent=4, sort_keys=True)
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
    run_id = 'a_run_id'
    _expected_output = None

    @property
    def expected_output(self):
        if self._expected_output is None:
            self._expected_output = json.load(
                open(os.path.join(self.test_data, 'expected_run_crawler_data.json'))
            )
        return self._expected_output

    def setUp(self):
        patched_lims_info = patch(ppath + 'RunCrawler.get_sample_information_from_lims')
        patched_data = patch(
            ppath + 'RunCrawler._run_sample_lane_to_barcode',
            return_value={
                'a_run_id_1_unknown': {'read_1_trimmed_bases': 184380158, 'read_2_trimmed_bases': 172552099},
                'a_run_id_2_unknown': {'read_1_trimmed_bases': 48149799, 'read_2_trimmed_bases': 48818739},
                'a_run_id_1_TCCGGAGA': {'read_1_trimmed_bases': 1088149481, 'read_2_trimmed_bases': 1034179505},
                'a_run_id_2_TCCGGAGA': {'read_1_trimmed_bases': 398993728, 'read_2_trimmed_bases': 391621660},
                'a_run_id_1_ATTACTCG': {'read_1_trimmed_bases': 714309214, 'read_2_trimmed_bases': 684692293},
                'a_run_id_2_ATTACTCG': {'read_1_trimmed_bases': 284712861, 'read_2_trimmed_bases': 282625840}}
        )
        run_element1 = {
            c.ELEMENT_PROJECT_ID: '10015AT',
            c.ELEMENT_SAMPLE_INTERNAL_ID: '10015AT0001',
            c.ELEMENT_LIBRARY_INTERNAL_ID: 'LP6002014-DTP_A01',
            c.ELEMENT_LANE: '1',
            c.ELEMENT_BARCODE: 'ATTACTCG'
        }
        run_element2 = {
            c.ELEMENT_PROJECT_ID: '10015AT',
            c.ELEMENT_SAMPLE_INTERNAL_ID: '10015AT0002',
            c.ELEMENT_LIBRARY_INTERNAL_ID: 'LP6002014-DTP_A02',
            c.ELEMENT_LANE: '1',
            c.ELEMENT_BARCODE: 'TCCGGAGA'
        }
        run_element3 = {
            c.ELEMENT_PROJECT_ID: '10015AT',
            c.ELEMENT_SAMPLE_INTERNAL_ID: '10015AT0001',
            c.ELEMENT_LIBRARY_INTERNAL_ID: 'LP6002014-DTP_A01',
            c.ELEMENT_LANE: '2',
            c.ELEMENT_BARCODE: 'ATTACTCG'
        }
        run_element4 = {
            c.ELEMENT_PROJECT_ID: '10015AT',
            c.ELEMENT_SAMPLE_INTERNAL_ID: '10015AT0002',
            c.ELEMENT_LIBRARY_INTERNAL_ID: 'LP6002014-DTP_A02',
            c.ELEMENT_LANE: '2',
            c.ELEMENT_BARCODE: 'TCCGGAGA'
        }
        dataset = NamedMock(
            real_name='a_run_id',
            run_elements=[run_element1, run_element2, run_element3, run_element4],
            has_barcodes=True
        )
        with patched_lims_info, patched_data:
            self.crawler = report_generation.RunCrawler(
                dataset,
                os.path.join(self.test_data, 'test_run_dir')
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

    def test_run_sample_lane_to_barcode(self):
        input_data = {(self.run_id, '10015AT0001', '1'): {'read_1_trimmed_bases': 714309214, 'read_2_trimmed_bases': 684692293}}
        test = self.crawler._run_sample_lane_to_barcode(input_data)
        assert test == {self.run_id + '_1_ATTACTCG': {'read_1_trimmed_bases': 714309214, 'read_2_trimmed_bases': 684692293}}


class TestSampleCrawler(TestCrawler):
    def setUp(self):
        self.expected_output = json.load(open(os.path.join(self.test_data, 'expected_sample_crawler_data.json')))
        patched_sample_info = patch(
            ppath + 'SampleCrawler.get_sample_information_from_lims',
            return_value={'user_sample_id': 'test_sample', 'provided_gender': 'female', 'species_name': 'Homo sapiens'}
        )
        patched_user_sample_id = patch(ppath + 'sample_crawler.clarity.get_user_sample_name', return_value='test_sample')
        output_cfg = OutputFileConfiguration('bcbio')
        with patched_sample_info, patched_user_sample_id:
            self.crawler = report_generation.SampleCrawler(
                'test_sample', 'test_project', self.test_data, output_cfg, post_pipeline=True
            )

    def test_sample(self):
        self.compare_jsons(self.crawler.sample, self.expected_output)
