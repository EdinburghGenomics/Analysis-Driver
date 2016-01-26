from tests.test_analysisdriver import TestAnalysisDriver
from analysis_driver.quality_control import ContaminationCheck
from analysis_driver.exceptions import AnalysisDriverError
from analysis_driver.config import default as cfg


class TestContaminationCheck(TestAnalysisDriver):

    def test_kontaminant_command1(self):
        fastq_files = ['fastqFile1', 'fastqFile2']
        sample_id = 'testSample'
        referenceSp = 'testSpecies'
        c = ContaminationCheck(fastq_files=fastq_files, sample_id=sample_id)
        kontaminantCmd = c._kontaminant_command(referenceSp)
        assert kontaminantCmd == "This is my paired end kontaminant command for testSpecies"

        fastq_files = ['fastqFile1']
        c = ContaminationCheck(fastq_files=fastq_files, sample_id=sample_id)
        kontaminantCmd = c._kontaminant_command(referenceSp)
        assert kontaminantCmd == "This is my single end kontaminant command for testSpecies"

        fastq_files = ['fastqFile1', 'fastqFile2', 'fastqFile3']
        c = ContaminationCheck(fastq_files=fastq_files, sample_id=sample_id)
        self.assertRaises(ValueError, c._kontaminant_command,referenceSp)


    def test_loop_through_references(self):
        fastq_files = ['fastqFile1', 'fastqFile2']
        sample_id = 'testSample'
        c = ContaminationCheck(fastq_files=fastq_files, sample_id=sample_id)
        referenceCommandDict = c._loop_through_references()
        assert referenceCommandDict == {'species1':'This is my paired end kontaminant command for species1',
                                        'species2':'This is my paired end kontaminant command for species2',
                                        'species3':'This is my paired end kontaminant command for species3'}












