__author__ = 'mwham'
from tests.test_analysisdriver import TestAnalysisDriver
from analysis_driver.reader import SampleSheet, RunInfo, transform_sample_sheet
from analysis_driver.reader.sample_sheet import SampleProject, Sample
import pytest
import os


samplesheet_with_barcode = """[Header],,,,,,,,
IEMFileVersion,4,,,,,,,
Investigator Name,Lucas Lefevre,,,,,,,
Experiment Name,,,,,,,,
Date,02/12/2015,,,,,,,
Workflow,GenerateFASTQ,,,,,,,
Application,HiSeq FASTQ Only,,,,,,,
Assay,TruSeq HT,,,,,,,
Description,HiSeqX run,,,,,,,
Chemistry,Default,,,,,,,
,,,,,,,,
[Reads],,,,,,,,
151,,,,,,,,
151,,,,,,,,
,,,,,,,,
[Settings],,,,,,,,
ReverseComplement,0,,,,,,,
Adapter,AGATCGGAAGAGCACACGTCTGAACTCCAGTCA,,,,,,,
AdapterRead2,AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT,,,,,,,
,,,,,,,,
[Data],
Lane,Sample_ID,Sample_Name,Sample_Plate,Sample_Well,Index2,Index,Sample_Project,GenomeFolder
1,10015AT0001,10015ATpool01,,,IL-TP-006,GCCAAT,10015AT,
2,10015AT0001,10015ATpool01,,,IL-TP-002,CGATGT,10015AT,
3,10015AT0001,10015ATpool01,,,IL-TP-007,CAGATC,10015AT,
4,10015AT0001,10015ATpool01,,,IL-TP-005,ACAGTG,10015AT,
5,10015AT0001,10015ATpool01,,,IL-TP-012,CTTGTA,10015AT,
6,10015AT0001,10015ATpool01,,,IL-TP-013,AGTCAA,10015AT,
7,10015AT0001,10015ATpool01,,,IL-TP-014,AGTTCC,10015AT,
8,10015AT0001,10015ATpool01,,,IL-TP-002,CGATGT,10015AT,
"""

samplesheet_without_barcode = """[Header],,,,,,,,
IEMFileVersion,4,,,,,,,
Investigator Name,Lucas Lefevre,,,,,,,
Experiment Name,,,,,,,,
Date,02/12/2015,,,,,,,
Workflow,GenerateFASTQ,,,,,,,
Application,HiSeq FASTQ Only,,,,,,,
Assay,TruSeq HT,,,,,,,
Description,HiSeqX run,,,,,,,
Chemistry,Default,,,,,,,
,,,,,,,,
[Reads],,,,,,,,
151,,,,,,,,
151,,,,,,,,
,,,,,,,,
[Settings],,,,,,,,
ReverseComplement,0,,,,,,,
Adapter,AGATCGGAAGAGCACACGTCTGAACTCCAGTCA,,,,,,,
AdapterRead2,AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT,,,,,,,
,,,,,,,,
[Data],
Lane,Sample_ID,Sample_Name,Sample_Plate,Sample_Well,Index2,Index,Sample_Project,GenomeFolder
1,10015AT0001,10015ATpool01,,,,,10015AT,
2,10015AT0001,10015ATpool01,,,,,10015AT,
3,10015AT0001,10015ATpool01,,,,,10015AT,
4,10015AT0001,10015ATpool01,,,,,10015AT,
5,10015AT0001,10015ATpool01,,,,,10015AT,
6,10015AT0001,10015ATpool01,,,,,10015AT,
7,10015AT0001,10015ATpool01,,,,,10015AT,
8,10015AT0001,10015ATpool01,,,,,10015AT,
"""
class TestSampleSheetFunctions(TestAnalysisDriver):
    new_SP_file = os.path.join(TestAnalysisDriver.assets_path, 'SampleSheet_analysis_driver.csv')

    def setUp(self):
        if os.path.exists(self.new_SP_file):
            os.unlink(self.new_SP_file)

    def tearDown(self):
        if os.path.exists(self.new_SP_file):
            os.unlink(self.new_SP_file)

    def test_transform_sample_sheet(self):
        assert not os.path.exists(self.new_SP_file)
        transform_sample_sheet(self.assets_path)
        assert os.path.exists(self.new_SP_file)
        with open(self.new_SP_file) as open_file:
            assert open_file.read() == samplesheet_with_barcode
        os.unlink(self.new_SP_file)

        transform_sample_sheet(self.assets_path, remove_barcode=True)
        assert os.path.exists(self.new_SP_file)
        with open(self.new_SP_file) as open_file:
            assert open_file.read() == samplesheet_without_barcode
        os.unlink(self.new_SP_file)


class TestSampleSheet(TestAnalysisDriver):
    def setUp(self):
        transform_sample_sheet(self.assets_path)
        self.sample_sheet = SampleSheet(self.sample_sheet_path)
        self.barcoded_samplesheet = SampleSheet(self.barcoded_samplesheet_path)
        self.barcodeless_samplesheet = SampleSheet(self.barcodeless_samplesheet_path)
        self.barcoded_run_info = RunInfo(self.barcoded_run_info_path)
        self.barcodeless_run_info = RunInfo(self.barcodeless_run_info_path)
        self.run_info = RunInfo(self.assets_path)
        self.samples = []
        for name, p in self.sample_sheet.sample_projects.items():
            for name2, i in p.sample_ids.items():
                for sample in i.samples:
                    self.samples.append(sample)
        self.sample_ids = self.sample_sheet.sample_projects['10015AT'].sample_ids

    def test_init(self):
        expected_lane = 1
        for sample in self.samples:
            assert sample.extra_data['Index2'] in [
                'IL-TP-002', 'IL-TP-005', 'IL-TP-006', 'IL-TP-007', 'IL-TP-012', 'IL-TP-013', 'IL-TP-014'
            ]
            assert sample.sample_project == '10015AT'
            assert sample.sample_id == '10015AT0001'
            assert sample.lane == str(expected_lane)
            assert sample.extra_data['Sample_Name'] == '10015ATpool01'
            assert sample.barcode == sample.barcode.upper()

            expected_lane += 1

    def test_check_barcodes(self):
        assert self.sample_sheet.check_barcodes() == 6

    def test_generate_mask(self):
        assert self.barcoded_samplesheet.generate_mask(self.barcoded_run_info.mask) == 'y150n,i8,y150n'
        assert self.barcodeless_samplesheet.generate_mask(self.barcodeless_run_info.mask) == 'y150n,y150n'

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
