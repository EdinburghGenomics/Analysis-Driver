from os.path import join
from tests.test_analysisdriver import TestAnalysisDriver
from analysis_driver.reader import SampleSheet, RunInfo, transform_sample_sheet
from analysis_driver.reader.sample_sheet import ProjectID, Line
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
    with open(new_samplesheet) as f:
        assert samples_without_barcode in f.read()

    transform_sample_sheet(samplesheet_path)
    with open(new_samplesheet) as f:
        assert samples_with_barcode in f.read()


def test_transform_sample_sheet_seqlab2():
    samplesheet_path = join(TestAnalysisDriver.assets_path, 'test_runs', 'clarity4_run')
    new_samplesheet = join(samplesheet_path, 'SampleSheet_analysis_driver.csv')

    transform_sample_sheet(samplesheet_path, clarity4=True)
    with open(new_samplesheet) as f:
        obs = f.read()
        assert samples_with_barcode in obs


class TestBarcodedSampleSheet(TestAnalysisDriver):
    def setUp(self):
        run_dir = join(self.assets_path, 'test_runs', 'barcoded_run')
        transform_sample_sheet(run_dir)

        self.sample_sheet = SampleSheet(join(run_dir, 'SampleSheet_analysis_driver.csv'))
        self.run_info = RunInfo(run_dir)
        self.samples = []
        for name, p in self.sample_sheet.projects.items():
            for name2, i in p.sample_ids.items():
                for sample in i.lines:
                    self.samples.append(sample)
        self.sample_ids = self.sample_sheet.projects['10015AT'].sample_ids

    def test_init(self):
        for sample in self.samples:
            assert sample.extra_data['Index2'] in [
                'IL-TP-002', 'IL-TP-005', 'IL-TP-006', 'IL-TP-007', 'IL-TP-012', 'IL-TP-013', 'IL-TP-014'
            ]
            assert sample.project_id == '10015AT'
            assert sample.sample_id == '10015AT0001'
            assert sample.lanes == ['1', '2', '3', '4', '5', '6', '7', '8']
            assert sample.sample_name == '10015ATpool01'
            assert sample.barcode == sample.barcode.upper()

    def test_check_barcodes(self):
        assert self.sample_sheet.barcode_len == 8

    def test_validate(self):
        assert self.sample_sheet.validate(self.run_info.reads) is True


class TestBarcodelessSampleSheet(TestAnalysisDriver):
    def setUp(self):
        run_dir = join(self.assets_path, 'test_runs', 'barcodeless_run')
        transform_sample_sheet(run_dir, remove_barcode=True)

        self.sample_sheet = SampleSheet(join(run_dir, 'SampleSheet_analysis_driver.csv'))
        self.run_info = RunInfo(run_dir)
        self.samples = []
        for name, p in self.sample_sheet.projects.items():
            for name2, i in p.sample_ids.items():
                for sample in i.lines:
                    self.samples.append(sample)
        self.sample_ids = self.sample_sheet.projects['10015AT'].sample_ids

    def test_init(self):
        expected_lane = 1
        for sample in self.samples:
            assert sample.project_id == '10015AT'
            assert sample.sample_id == '10015AT0001'
            assert sample.lanes == [str(expected_lane)]
            assert sample.sample_name == '10015ATpool01'
            assert sample.barcode == sample.barcode.upper()

            expected_lane += 1

    def test_check_barcodes(self):
        assert self.sample_sheet.barcode_len == 0
        assert self.sample_sheet.has_barcodes is False

    def test_check_one_barcode_per_lane(self):
        assert self.sample_sheet._validate_one_sample_per_lane() is None  # no error raised

    def test_validate(self):
        assert self.sample_sheet.validate(self.run_info.reads) is True


class TestSampleProject(TestAnalysisDriver):
    def setUp(self):
        self.test_line = Line(
            {
                'Sample_Project': 'test_sp',
                'Lane': '1337',
                'Sample_ID': 'test_id',
                'Sample_Name': 'test_name',
                'Index': 'ATGCAT'
            }
        )
        self.sample_project = ProjectID('test_sp')
        self.sample_project.get_sample_id(self.test_line.sample_id).add_line(self.test_line)

    def test_init(self):
        assert self.sample_project.name == 'test_sp'
        assert self.sample_project.sample_ids['test_id'].lines[0] == self.test_line

    def test_add_sample(self):
        new_line = Line(
            {
                'Sample_Project': 'test_sp',
                'Lane': '1338',
                'Sample_ID': 'test_id',
                'Sample_Name': 'test_name',
                'Index': 'ATGCAG'
            }
        )
        self.sample_project.get_sample_id('test_id').add_line(new_line)

        assert new_line != self.test_line
        assert new_line in self.sample_project.get_sample_id('test_id').lines

    def test_add_faulty_sample(self):
        with pytest.raises(AssertionError) as e:
            new_line = Line(
                {
                    'Sample_Project': 'another_test_sp',
                    'Lane': '1338',
                    'Sample_ID': 'another_test_id',
                    'Sample_Name': 'test_name',
                    'Index': 'ATGCAG'
                }
            )
            self.sample_project.get_sample_id('test_id').add_line(new_line)
        assert 'Adding invalid sample project to test_id: another_test_sp' == str(e.value)
