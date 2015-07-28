__author__ = 'mwham'
from tests.test_analysisdriver import TestAnalysisDriver
from analysis_driver.reader.sample_sheet import SampleSheet, SampleProject, Sample
import pytest


class TestSampleSheet(TestAnalysisDriver):
    def setUp(self):
        self.sample_sheet = SampleSheet(self.assets_path)
        self.samples = self.sample_sheet.sample_projects['10015AT'].samples

    def test_init(self):
        expected_lane = 1
        for sample in self.samples:
            assert sample.index_id in [
                'IL-TP-002', 'IL-TP-005', 'IL-TP-006', 'IL-TP-007', 'IL-TP-012', 'IL-TP-013', 'IL-TP-014'
            ]
            assert sample.sample_project == '10015AT'
            assert sample.id == '10015TA0001L05'
            assert sample.lane == str(expected_lane)
            assert sample.name == '10015ATpool01'
            assert sample.barcode == sample.barcode.upper()

            expected_lane += 1

    def test_check_barcodes(self):
        assert self.sample_sheet.check_barcodes() == 6


class TestSampleProject(TestAnalysisDriver):
    def setUp(self):
        self.test_sample = Sample(
            sample_project='test_sp',
            lane='1337',
            id='test_id',
            name='test_name',
            barcode='ATGCAT'
        )
        self.sample_project = SampleProject('test_sp', self.test_sample)

    def test_init(self):
        assert self.sample_project.name == 'test_sp'
        assert self.sample_project.samples[0] == self.test_sample

    def test_add_sample(self):
        new_sample = Sample(
            sample_project='test_sp',
            lane='1338',
            id='another_test_id',
            name='test_name',
            barcode='ATGCAG'
        )
        self.sample_project.add_sample(new_sample)

        assert self.sample_project.samples[1] == new_sample
        assert new_sample != self.test_sample

    def test_add_faulty_sample(self):
        with pytest.raises(AssertionError) as e:
            new_sample = Sample(
                sample_project='another_test_sp',
                lane='1338',
                id='another_test_id',
                name='another_test_name',
                barcode='ATGCAG'
            )
            self.sample_project.add_sample(new_sample)
        assert 'Adding invalid sample project to test_sp: another_test_sp' == str(e.value)
