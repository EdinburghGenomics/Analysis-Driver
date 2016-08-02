from tests.test_quality_control.qc_tester import QCTester
from analysis_driver.quality_control import ContaminationCheck, VerifyBamId
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

    @patch('analysis_driver.dataset.rest_communication')
    @patch('egcg_core.executor.execute')
    def test_run_fastqscreen(self, mocked_execute, mocked_rest):
        instance = mocked_execute.return_value
        instance.join.return_value = 0
        run_fastqscreen_single_end = self.c_se._run_fastqscreen()

        assert run_fastqscreen_single_end == ['fastqFile1_screen.txt']
        mocked_execute.assert_called_once_with(
            'path/to/fastqscreen --aligner bowtie2 fastqFile1.fastq '
            '--conf path/to/fastqscreen/conf --force',
            working_dir='testSample',
            mem=10,
            cpus=2,
            job_name='fastqscreen'
        )
        mocked_execute.reset_mock()
        run_fastqscreen_paired_end = self.c_pe._run_fastqscreen()
        assert run_fastqscreen_paired_end == ['fastqFile1_screen.txt', 'fastqFile2_screen.txt']
        mocked_execute.assert_called_once_with(
            'path/to/fastqscreen --aligner bowtie2 '
            'fastqFile1.fastq fastqFile2.fastq '
            '--conf path/to/fastqscreen/conf --force',
            working_dir='testSample',
            mem=10,
            cpus=2,
            job_name='fastqscreen'
        )


class TestVerifyBamId(QCTester):
    def setUp(self):
        super().setUp()
        self.working_dir = 'testSample'
        self.vbi = VerifyBamId(self.dataset, self.working_dir, 'test_bam_file.bam')

    @patch('analysis_driver.dataset.rest_communication')
    @patch('egcg_core.executor.execute')
    def test_contamination_check(self, mocked_execute, mocked_rest):
        self.vbi._contamination_check()
        assert mocked_execute.call_count == 3
        mocked_execute.assert_any_call(
            'path/to/samtools view -b test_bam_file.bam chr22 > testSample/test_sample_chr22.bam',
            mem=2,
            job_name='filter_bam22',
            working_dir='testSample',
            cpus=1,
            log_commands=False
        )
        mocked_execute.assert_any_call(
                'path/to/samtools index testSample/test_sample_chr22.bam',
                mem=2,
                working_dir='testSample',
                cpus=1,
                job_name='index_bam22'
        )
        mocked_execute.assert_any_call(
            'path/to/verifybamid --bam testSample/test_sample_chr22.bam --vcf path/to/population_vcf --out testSample/test_sample-chr22-vbi',
            working_dir='testSample',
            cpus=1,
            job_name='verify_bam_id',
            mem=4
        )


    @patch('egcg_core.executor.execute')
    def test_filter_bam(self, mocked_execute):
        self.vbi._filter_bam()
        mocked_execute.assert_called_once_with(
            'path/to/samtools view -b test_bam_file.bam chr22 > testSample/test_sample_chr22.bam',
            job_name='filter_bam22',
            mem=2,
            working_dir='testSample',
            cpus=1,
            log_commands=False
        )

    @patch('egcg_core.executor.execute')
    def test_index_filtered_bam(self, mocked_execute):
        self.vbi.filtered_bam = 'test_filtered_bam.bam'
        self.vbi._index_filtered_bam()
        mocked_execute.assert_called_once_with(
            'path/to/samtools index test_filtered_bam.bam',
            mem=2,
            job_name='index_bam22',
            working_dir='testSample',
            cpus=1
        )

    @patch('egcg_core.executor.execute')
    def test_verify_bam_id(self, mocked_execute):
        self.vbi.filtered_bam = 'test_filtered_bam.bam'
        self.vbi._verify_bam_id()
        mocked_execute.assert_called_once_with(
            'path/to/verifybamid --bam test_filtered_bam.bam --vcf path/to/population_vcf --out testSample/test_sample-chr22-vbi',
            job_name='verify_bam_id',
            mem=4,
            cpus=1,
            working_dir='testSample'
        )
