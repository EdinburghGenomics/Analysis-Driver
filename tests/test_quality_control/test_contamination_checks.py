from tests.test_analysisdriver import TestAnalysisDriver
from analysis_driver.quality_control import ContaminationCheck
from analysis_driver.config import default as cfg
from unittest.mock import patch


class TestContaminationCheck(TestAnalysisDriver):

    def test_createOutputDirectories(self):
        fastq_files = ['fastqFile1', 'fastqFile2']
        sample_id = 'testSample'
        species = 'testSpecies'
        c = ContaminationCheck(fastq_files=fastq_files, sample_id=sample_id)
        createOutputDirectoriesCommand = c.createOutputDirectories(species)
        assert createOutputDirectoriesCommand =='path/to/kontaminant/working/dir/testSample/testSpecies'




    def test_kontaminant_command1(self):
        fastq_files = ['fastqFile1', 'fastqFile2']
        sample_id = 'testSample'
        referenceSp = 'testSpecies'
        outputDir = 'testOutputDir'
        c = ContaminationCheck(fastq_files=fastq_files, sample_id=sample_id)
        kontaminantCmd = c._kontaminant_command(referenceSp, outputDir)
        assert kontaminantCmd == "path/to/kontaminant/bin/kmer_filter_31 --reference path/to/kontaminant/reference/databases/testSpecies --read1 fastqFile1 --read2 fastqFile2 --output_prefix filtered_K21_Thr1 --output_folder testOutputDir --threashold 1"

        fastq_files = ['fastqFile1']
        c = ContaminationCheck(fastq_files=fastq_files, sample_id=sample_id)
        kontaminantCmd = c._kontaminant_command(referenceSp, outputDir)
        assert kontaminantCmd == "path/to/kontaminant/bin/kmer_filter_31 --reference path/to/kontaminant/reference/databases/testSpecies --read1 fastqFile1 --output_prefix filtered_K21_Thr1 --output_folder testOutputDir --threashold 1"

        fastq_files = ['fastqFile1', 'fastqFile2', 'fastqFile3']
        c = ContaminationCheck(fastq_files=fastq_files, sample_id=sample_id)
        self.assertRaises(ValueError, c._kontaminant_command,referenceSp,outputDir)


    def test_loop_through_references(self):
        fastq_files = ['fastqFile1', 'fastqFile2']
        sample_id = 'testSample'
        c = ContaminationCheck(fastq_files=fastq_files, sample_id=sample_id)
        referenceCommandDict = c._loop_through_references()
        print(referenceCommandDict)
        assert referenceCommandDict == {'species1':'path/to/kontaminant/bin/kmer_filter_31 --reference path/to/kontaminant/reference/databases/species1 --read1 fastqFile1 --read2 fastqFile2 --output_prefix filtered_K21_Thr1 --output_folder path/to/kontaminant/working/dir/testSample/species1 --threashold 1',
                                        'species2':'path/to/kontaminant/bin/kmer_filter_31 --reference path/to/kontaminant/reference/databases/species2 --read1 fastqFile1 --read2 fastqFile2 --output_prefix filtered_K21_Thr1 --output_folder path/to/kontaminant/working/dir/testSample/species2 --threashold 1',
                                        'species3':'path/to/kontaminant/bin/kmer_filter_31 --reference path/to/kontaminant/reference/databases/species3 --read1 fastqFile1 --read2 fastqFile2 --output_prefix filtered_K21_Thr1 --output_folder path/to/kontaminant/working/dir/testSample/species3 --threashold 1'}

    @patch('analysis_driver.executor.execute', autospec=True)
    def test_runKontaminant(self, mocked_execute):
        fastq_files = ['fastqFile1', 'fastqFile2']
        sample_id = 'testSample'
        c = ContaminationCheck(fastq_files=fastq_files, sample_id=sample_id)
        instance = mocked_execute.return_value
        instance.join.return_value = 0
        runContaminantCommand = c._runKontaminant()
        assert runContaminantCommand == {'species1':"path/to/kontaminant/working/dir/testSample/species1/filtered_K21_Thr1species1.FASTQ",
                                         'species2':"path/to/kontaminant/working/dir/testSample/species2/filtered_K21_Thr1species2.FASTQ",
                                         'species3':"path/to/kontaminant/working/dir/testSample/species3/filtered_K21_Thr1species3.FASTQ"}











