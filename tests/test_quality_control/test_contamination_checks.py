from tests.test_quality_control.qc_tester import QCTester
from analysis_driver.quality_control import ContaminationCheck, VerifyBamID
from unittest.mock import patch


class TestContaminationCheck(QCTester):
    def setUp(self):
        super().setUp()
        self.c_se = ContaminationCheck(dataset=self.dataset, fastq_files=['sample_r1.fastq.gz'])
        self.c_pe = ContaminationCheck(dataset=self.dataset, fastq_files=['sample_r1.fastq.gz', 'sample_r2.fastq.gz'])

    def test_fastqscreen_command(self):
        assert self.c_se._fastqscreen_command() == (
            'path/to/fastqscreen --aligner bowtie2 sample_r1.fastq.gz --conf path/to/fastqscreen/conf --force'
        )
        assert self.c_pe._fastqscreen_command() == ('path/to/fastqscreen --aligner bowtie2 sample_r1.fastq.gz'
                                                    ' sample_r2.fastq.gz --conf path/to/fastqscreen/conf --force')

    def test_get_expected_outfiles(self):
        assert self.c_se.fastqscreen_expected_outfiles == ['sample_r1_screen.txt']
        assert self.c_pe.fastqscreen_expected_outfiles == ['sample_r1_screen.txt', 'sample_r2_screen.txt']

    @patch('egcg_core.executor.execute')
    def test_run_fastqscreen(self, mocked_execute):
        self.c_se._run()
        mocked_execute.assert_called_once_with(
            'path/to/fastqscreen --aligner bowtie2 sample_r1.fastq.gz '
            '--conf path/to/fastqscreen/conf --force',
            working_dir='path/to/jobs/test_sample',
            mem=10,
            cpus=2,
            job_name='fastqscreen'
        )
        mocked_execute.reset_mock()
        self.c_pe._run()
        mocked_execute.assert_called_once_with(
            'path/to/fastqscreen --aligner bowtie2 '
            'sample_r1.fastq.gz sample_r2.fastq.gz '
            '--conf path/to/fastqscreen/conf --force',
            working_dir='path/to/jobs/test_sample',
            mem=10,
            cpus=2,
            job_name='fastqscreen'
        )


class TestVerifyBamId(QCTester):
    def setUp(self):
        super().setUp()
        self.vbi = VerifyBamID(dataset=self.dataset, bam_file='test_bam_file.bam')

    @patch('egcg_core.executor.execute')
    def test_contamination_check(self, mocked_execute):
        self.vbi._run()
        mocked_execute.assert_any_call(
            'path/to/samtools view -b test_bam_file.bam chr22 > path/to/jobs/test_sample/test_sample_chr22.bam',
            mem=2,
            job_name='filter_bam22',
            working_dir='path/to/jobs/test_sample',
            cpus=1,
            log_commands=False
        )
        mocked_execute.assert_any_call(
            'path/to/samtools index path/to/jobs/test_sample/test_sample_chr22.bam',
            mem=2,
            working_dir='path/to/jobs/test_sample',
            cpus=1,
            job_name='index_bam22'
        )
        mocked_execute.assert_any_call(
            'path/to/verifybamid --bam path/to/jobs/test_sample/test_sample_chr22.bam --vcf path/to/population_vcf --out path/to/jobs/test_sample/test_sample-chr22-vbi',
            working_dir='path/to/jobs/test_sample',
            cpus=1,
            job_name='verify_bam_id',
            mem=4
        )

    @patch('egcg_core.executor.execute')
    def test_filter_bam(self, mocked_execute):
        self.vbi._filter_bam()
        mocked_execute.assert_called_once_with(
            'path/to/samtools view -b test_bam_file.bam chr22 > path/to/jobs/test_sample/test_sample_chr22.bam',
            job_name='filter_bam22',
            mem=2,
            working_dir='path/to/jobs/test_sample',
            cpus=1,
            log_commands=False
        )

    @patch('egcg_core.executor.execute')
    def test_index_filtered_bam(self, mocked_execute):
        self.vbi._index_filtered_bam()
        mocked_execute.assert_called_once_with(
            'path/to/samtools index path/to/jobs/test_sample/test_sample_chr22.bam',
            mem=2,
            job_name='index_bam22',
            working_dir='path/to/jobs/test_sample',
            cpus=1
        )

    @patch('egcg_core.executor.execute')
    def test_verify_bam_id(self, mocked_execute):
        self.vbi._verify_bam_id()
        mocked_execute.assert_called_once_with(
            'path/to/verifybamid --bam path/to/jobs/test_sample/test_sample_chr22.bam --vcf path/to/population_vcf --out path/to/jobs/test_sample/test_sample-chr22-vbi',
            job_name='verify_bam_id',
            mem=4,
            cpus=1,
            working_dir='path/to/jobs/test_sample'
        )
