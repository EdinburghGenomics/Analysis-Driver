from tests.test_analysisdriver import TestAnalysisDriver
from analysis_driver.quality_control import ContaminationCheck
from unittest.mock import patch


class TestContaminationCheck(TestAnalysisDriver):

    def test_fastqscreen_command(self):
        fastq_files = ['fastqFile1.fastq', 'fastqFile2.fastq']
        working_dir = 'testSample'
        c = ContaminationCheck(fastq_files=fastq_files, working_dir=working_dir)
        my_fastq_screen_command = c._fastqscreen_command()
        assert my_fastq_screen_command == ["path/to/fastqscreen " \
                                          "--aligner bowtie2 fastqFile1.fastq fastqFile2.fastq " \
                                          "--conf path/to/fastqscreen/conf --force"]

    def test_get_expected_outfiles(self):
        fastq_files = ['fastqFile1.fastq', 'fastqFile2.fastq']
        working_dir = 'testSample'
        c = ContaminationCheck(fastq_files=fastq_files, working_dir=working_dir)
        fastqscreen_expected_outfiles = c._get_expected_outfiles()
        assert fastqscreen_expected_outfiles == ['fastqFile1_screen.txt', 'fastqFile2_screen.txt']

    @patch('analysis_driver.executor.execute', autospec=True)
    def test_run_fastqscreen(self, mocked_execute):
        fastq_files = ['fastqFile1.fastq', 'fastqFile2.fastq']
        working_dir = 'testSample'
        c = ContaminationCheck(fastq_files=fastq_files, working_dir=working_dir)
        instance = mocked_execute.return_value
        instance.join.return_value = 0
        run_fastqscreen = c._run_fastqscreen()
        assert run_fastqscreen == ['fastqFile1_screen.txt', 'fastqFile2_screen.txt']

    @patch('analysis_driver.executor.execute', autospec=True)
    def test_run_fastqscreen(self, mocked_execute):
        fastq_files = ['fastqFile1.fastq', 'fastqFile2.fastq']
        working_dir = 'testSample'
        c = ContaminationCheck(fastq_files=fastq_files, working_dir=working_dir)
        instance = mocked_execute.return_value
        instance.join.return_value = 0
        run_fastqscreen = c._run_fastqscreen()
        mocked_execute.assert_called_once_with(['path/to/fastqscreen '
                                                '--aligner bowtie2 fastqFile1.fastq fastqFile2.fastq '
                                                '--conf path/to/fastqscreen/conf --force'],
                                               working_dir='testSample',
                                               mem=10,
                                               cpus=2,
                                               job_name='fastqscreen')



