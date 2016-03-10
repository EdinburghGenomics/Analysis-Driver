from tests.test_analysisdriver import TestAnalysisDriver
from analysis_driver.quality_control import ContaminationCheck
from unittest.mock import patch


class TestContaminationCheck(TestAnalysisDriver):

    def test_fastqscreen_command_pairedend(self):
        fastq_files = ['fastqFile1.fastq', 'fastqFile2.fastq']
        working_dir = 'testSample'
        c = ContaminationCheck(fastq_files=fastq_files, working_dir=working_dir)
        my_fastq_screen_command = c._fastqscreen_command()
        assert my_fastq_screen_command == ("path/to/fastqscreen " \
                                          "--aligner bowtie2 fastqFile1.fastq fastqFile2.fastq " \
                                          "--conf path/to/fastqscreen/conf --force")

    def test_fastqscreen_command_singleend(self):
        fastq_files = ['fastqFile1.fastq']
        working_dir = 'testSample'
        c = ContaminationCheck(fastq_files=fastq_files, working_dir=working_dir)
        my_fastq_screen_command = c._fastqscreen_command()
        assert my_fastq_screen_command == ("path/to/fastqscreen " \
                                          "--aligner bowtie2 fastqFile1.fastq " \
                                          "--conf path/to/fastqscreen/conf --force")

    def test_get_expected_outfiles(self):
        fastq_files = ['fastqFile1.fastq', 'fastqFile2.fastq']
        working_dir = 'testSample'
        c = ContaminationCheck(fastq_files=fastq_files, working_dir=working_dir)
        fastqscreen_expected_outfiles = c._get_expected_outfiles()
        assert fastqscreen_expected_outfiles == ['fastqFile1_screen.txt', 'fastqFile2_screen.txt']

    @patch('analysis_driver.executor.execute')
    def test_run_fastqscreen(self, mocked_execute):
        fastq_files_single_end = ['fastqFile1.fastq']
        fastq_files_paired_end = ['fastqFile1.fastq', 'fastqFile2.fastq']
        working_dir = 'testSample'
        c_single_end = ContaminationCheck(fastq_files=fastq_files_single_end, working_dir=working_dir)
        c_paired_end = ContaminationCheck(fastq_files=fastq_files_paired_end, working_dir=working_dir)
        instance = mocked_execute.return_value
        instance.join.return_value = 0
        run_fastqscreen_single_end = c_single_end._run_fastqscreen()

        assert run_fastqscreen_single_end == ['fastqFile1_screen.txt']
        mocked_execute.assert_called_once_with(['path/to/fastqscreen '
                                                '--aligner bowtie2 fastqFile1.fastq '
                                                '--conf path/to/fastqscreen/conf --force'],
                                               working_dir='testSample',
                                               mem=10,
                                               cpus=2,
                                               job_name='fastqscreen')
        mocked_execute.reset_mock()
        run_fastqscreen_paired_end = c_paired_end._run_fastqscreen()
        assert run_fastqscreen_paired_end == ['fastqFile1_screen.txt', 'fastqFile2_screen.txt']
        mocked_execute.assert_called_once_with(['path/to/fastqscreen '
                                                '--aligner bowtie2 fastqFile1.fastq fastqFile2.fastq '
                                                '--conf path/to/fastqscreen/conf --force'],
                                               working_dir='testSample',
                                               mem=10,
                                               cpus=2,
                                               job_name='fastqscreen')


