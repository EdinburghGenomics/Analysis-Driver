from tests.test_quality_control.qc_tester import QCTester
from analysis_driver.quality_control import ContaminationCheck
from unittest.mock import patch


class TestContaminationCheck(QCTester):
    def setUp(self):
        super().setUp()
        self.fastq_files_single_end = ['fastqFile1.fastq']
        self.fastq_files_paired_end = ['fastqFile1.fastq', 'fastqFile2.fastq']
        self.working_dir = 'testSample'
        self.c_se = ContaminationCheck(self.dataset, self.working_dir, self.fastq_files_single_end)
        self.c_pe = ContaminationCheck(self.dataset, self.working_dir, self.fastq_files_paired_end)

    def test_fastqscreen_command(self):
        my_fastq_screen_command_se = self.c_se._fastqscreen_command()
        my_fastq_screen_command_pe = self.c_pe._fastqscreen_command()
        assert my_fastq_screen_command_se == ('path/to/fastqscreen '
                                              '--aligner bowtie2 fastqFile1.fastq '
                                              '--conf path/to/fastqscreen/conf --force')
        assert my_fastq_screen_command_pe == ('path/to/fastqscreen '
                                              '--aligner bowtie2 fastqFile1.fastq fastqFile2.fastq '
                                              '--conf path/to/fastqscreen/conf --force')

    def test_get_expected_outfiles(self):
        fastqscreen_expected_outfiles_se = self.c_se._get_expected_outfiles()
        fastqscreen_expected_outfiles_pe = self.c_pe._get_expected_outfiles()
        assert fastqscreen_expected_outfiles_se == ['fastqFile1_screen.txt']
        assert fastqscreen_expected_outfiles_pe == ['fastqFile1_screen.txt', 'fastqFile2_screen.txt']

    @patch('analysis_driver.dataset_scanner.rest_communication')
    @patch('analysis_driver.executor.execute')
    def test_run_fastqscreen(self, mocked_execute, mocked_rest):
        instance = mocked_execute.return_value
        instance.join.return_value = 0
        run_fastqscreen_single_end = self.c_se._run_fastqscreen()

        assert run_fastqscreen_single_end == ['fastqFile1_screen.txt']
        mocked_execute.assert_called_once_with(
            ['path/to/fastqscreen --aligner bowtie2 fastqFile1.fastq '
             '--conf path/to/fastqscreen/conf --force'],
            working_dir='testSample',
            mem=10,
            cpus=2,
            job_name='fastqscreen'
        )
        mocked_execute.reset_mock()
        run_fastqscreen_paired_end = self.c_pe._run_fastqscreen()
        assert run_fastqscreen_paired_end == ['fastqFile1_screen.txt', 'fastqFile2_screen.txt']
        mocked_execute.assert_called_once_with(
            ['path/to/fastqscreen --aligner bowtie2 '
             'fastqFile1.fastq fastqFile2.fastq '
             '--conf path/to/fastqscreen/conf --force'],
            working_dir='testSample',
            mem=10,
            cpus=2,
            job_name='fastqscreen'
        )
