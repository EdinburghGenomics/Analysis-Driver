from tests.test_analysisdriver import TestAnalysisDriver
from analysis_driver.quality_control import ContaminationCheck
from analysis_driver.config import default as cfg
from unittest.mock import patch


class TestContaminationCheck(TestAnalysisDriver):

    def test_fastqscreen_command(self):
        fastq_files = ['fastqFile1.fastq', 'fastqFile2.fastq']
        sample_id = 'testSample'
        c = ContaminationCheck(fastq_files=fastq_files, sample_id=sample_id)
        my_fastq_screen_command = c._fastqscreen_command()
        assert my_fastq_screen_command == ["path/to/fastq-screen/bin " \
                                          "--aligner bowtie2 fastqFile1.fastq fastqFile2.fastq " \
                                          "--conf path/to/fastq-screen/conf/file"]

    def test_get_expected_outfiles(self):
        fastq_files = ['fastqFile1.fastq', 'fastqFile2.fastq']
        sample_id = 'testSample'
        c = ContaminationCheck(fastq_files=fastq_files, sample_id=sample_id)
        fastqscreen_expected_outfiles = c._get_expected_outfiles()
        print(fastqscreen_expected_outfiles)
        assert fastqscreen_expected_outfiles == ['fastqFile1_screen.txt', 'fastqFile2_screen.txt']

    @patch('analysis_driver.executor.execute', autospec=True)
    def test_run_fastqscreen(self, mocked_execute):
        fastq_files = ['fastqFile1.fastq', 'fastqFile2.fastq']
        sample_id = 'testSample'
        c = ContaminationCheck(fastq_files=fastq_files, sample_id=sample_id)
        instance = mocked_execute.return_value
        instance.join.return_value = 0
        run_fastqscreen = c._run_fastqscreen()
        assert run_fastqscreen == ['fastqFile1_screen.txt', 'fastqFile2_screen.txt']





