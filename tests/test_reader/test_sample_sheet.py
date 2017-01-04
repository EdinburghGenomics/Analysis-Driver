from os.path import join, isfile
from tests.test_analysisdriver import TestAnalysisDriver
from analysis_driver.reader import SampleSheet, RunInfo, transform_sample_sheet
from analysis_driver.reader.sample_sheet import SampleProject, Sample
import pytest


samples_with_barcode = """
Lane,Sample_ID,Sample_Name,Sample_Plate,Sample_Well,Index2,Index,Sample_Project,GenomeFolder
1+2+3+4+5+6+7+8,10015AT0001,10015ATpool01,,,IL-TP-006,GCCAATAA,10015AT,
1+2+3+4+5+6+7+8,10015AT0001,10015ATpool01,,,IL-TP-002,CGATGTAA,10015AT,
1+2+3+4+5+6+7+8,10015AT0001,10015ATpool01,,,IL-TP-007,CAGATCAA,10015AT,
1+2+3+4+5+6+7+8,10015AT0001,10015ATpool01,,,IL-TP-005,ACAGTGAA,10015AT,
1+2+3+4+5+6+7+8,10015AT0001,10015ATpool01,,,IL-TP-012,CTTGTAAA,10015AT,
1+2+3+4+5+6+7+8,10015AT0001,10015ATpool01,,,IL-TP-013,AGTCAAAA,10015AT,
1+2+3+4+5+6+7+8,10015AT0001,10015ATpool01,,,IL-TP-014,AGTTCCAA,10015AT,
1+2+3+4+5+6+7+8,10015AT0001,10015ATpool01,,,IL-TP-002,CGATGTAA,10015AT,
"""

samples_without_barcode = """
Lane,Sample_ID,Sample_Name,Sample_Plate,Sample_Well,Index2,Index,Sample_Project,GenomeFolder
1+2+3+4+5+6+7+8,10015AT0001,10015ATpool01,,,,,10015AT,
1+2+3+4+5+6+7+8,10015AT0001,10015ATpool01,,,,,10015AT,
1+2+3+4+5+6+7+8,10015AT0001,10015ATpool01,,,,,10015AT,
1+2+3+4+5+6+7+8,10015AT0001,10015ATpool01,,,,,10015AT,
1+2+3+4+5+6+7+8,10015AT0001,10015ATpool01,,,,,10015AT,
1+2+3+4+5+6+7+8,10015AT0001,10015ATpool01,,,,,10015AT,
1+2+3+4+5+6+7+8,10015AT0001,10015ATpool01,,,,,10015AT,
1+2+3+4+5+6+7+8,10015AT0001,10015ATpool01,,,,,10015AT,
"""


def test_transform_sample_sheet():
    samplesheet_path = join(TestAnalysisDriver.assets_path, 'test_runs', 'barcoded_run')
    new_samplesheet = join(samplesheet_path, 'SampleSheet_analysis_driver.csv')

    transform_sample_sheet(samplesheet_path, remove_barcode=True)
    assert isfile(new_samplesheet)
    with open(new_samplesheet) as open_file:
        obs = open_file.read()
        assert samples_without_barcode in obs

    transform_sample_sheet(samplesheet_path)
    assert isfile(new_samplesheet)
    with open(new_samplesheet) as open_file:
        obs = open_file.read()
        assert samples_with_barcode in obs


class TestBarcodedSampleSheet(TestAnalysisDriver):
    def setUp(self):
        run_dir = join(self.assets_path, 'test_runs', 'barcoded_run')
        transform_sample_sheet(run_dir)

        self.sample_sheet = SampleSheet(join(run_dir, 'SampleSheet_analysis_driver.csv'))
        self.run_info = RunInfo(run_dir)
        self.samples = []
        for name, p in self.sample_sheet.sample_projects.items():
            for name2, i in p.sample_ids.items():
                for sample in i.samples:
                    self.samples.append(sample)
        self.sample_ids = self.sample_sheet.sample_projects['10015AT'].sample_ids

    def test_init(self):
        for sample in self.samples:
            assert sample.extra_data['Index2'] in [
                'IL-TP-002', 'IL-TP-005', 'IL-TP-006', 'IL-TP-007', 'IL-TP-012', 'IL-TP-013', 'IL-TP-014'
            ]
            assert sample.sample_project == '10015AT'
            assert sample.sample_id == '10015AT0001'
            assert sample.lane == '1+2+3+4+5+6+7+8'
            assert sample.extra_data['Sample_Name'] == '10015ATpool01'
            assert sample.barcode == sample.barcode.upper()

    def test_check_barcodes(self):
        assert self.sample_sheet.check_barcodes() == 8

    def test_generate_mask(self):
        assert self.sample_sheet.generate_mask(self.run_info.mask) == 'y150n,i8,y150n'

    def test_validate(self):
        assert self.sample_sheet.validate(self.run_info.mask) is True


class TestBarcodelessSampleSheet(TestAnalysisDriver):
    def setUp(self):
        run_dir = join(self.assets_path, 'test_runs', 'barcodeless_run')

        transform_sample_sheet(run_dir, remove_barcode=True)

        self.sample_sheet = SampleSheet(join(run_dir, 'SampleSheet_analysis_driver.csv'), has_barcode=False)
        self.run_info = RunInfo(run_dir)
        self.samples = []
        for name, p in self.sample_sheet.sample_projects.items():
            for name2, i in p.sample_ids.items():
                for sample in i.samples:
                    self.samples.append(sample)
        self.sample_ids = self.sample_sheet.sample_projects['10015AT'].sample_ids

    def test_init(self):
        expected_lane = 1
        for sample in self.samples:
            assert sample.sample_project == '10015AT'
            assert sample.sample_id == '10015AT0001'
            assert sample.lane == str(expected_lane)
            assert sample.extra_data['Sample_Name'] == '10015ATpool01'
            assert sample.barcode == sample.barcode.upper()

            expected_lane += 1

    def test_check_barcodes(self):
        assert self.sample_sheet.check_barcodes() == 0

    def test_generate_mask(self):
        assert self.sample_sheet.generate_mask(self.run_info.mask) == 'y150n,y150n'

    def test_check_one_barcode_per_lane(self):
        self.sample_sheet._validate_one_sample_per_lane()

    def test_validate(self):
        assert self.sample_sheet.validate(self.run_info.mask) is True


class TestSampleProject(TestAnalysisDriver):
    def setUp(self):
        self.test_sample = Sample(
            sample_project='test_sp',
            lane='1337',
            sample_id='test_id',
            sample_name='test_name',
            name='test_name',
            barcode='ATGCAT'
        )
        self.sample_project = SampleProject('test_sp')
        self.sample_project.get_sample_id(self.test_sample.sample_id).add_sample(self.test_sample)

    def test_init(self):
        assert self.sample_project.name == 'test_sp'
        assert self.sample_project.sample_ids['test_id'].samples[0] == self.test_sample

    def test_add_sample(self):
        new_sample = Sample(
            sample_project='test_sp',
            lane='1338',
            sample_id='test_id',
            sample_name='test_name',
            name='test_name',
            barcode='ATGCAG'
        )
        self.sample_project.get_sample_id('test_id').add_sample(new_sample)

        assert new_sample != self.test_sample
        assert new_sample in self.sample_project.get_sample_id('test_id').samples

    def test_add_faulty_sample(self):
        with pytest.raises(AssertionError) as e:
            new_sample = Sample(
                sample_project='another_test_sp',
                lane='1338',
                sample_id='another_test_id',
                sample_name='test_name',
                name='another_test_name',
                barcode='ATGCAG'
            )
            self.sample_project.get_sample_id('test_id').add_sample(new_sample)
        assert 'Adding invalid sample project to test_id: another_test_sp' == str(e.value)
