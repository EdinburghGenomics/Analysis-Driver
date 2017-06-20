from tests.test_quality_control.qc_tester import QCTester
from analysis_driver.quality_control import ContaminationCheck, VerifyBamID
from unittest.mock import patch

ppath = 'analysis_driver.quality_control.contamination_checks.'


class TestContaminationCheck(QCTester):
    def setUp(self):
        super().setUp()
        self.c = ContaminationCheck(dataset=self.dataset, fq_pattern='r?.fastq.gz')

    def test_fastqscreen_command(self):
        with patch(ppath + 'util.find_files', return_value=['r1.fastq.gz', 'r2.fastq.gz']):
            assert self.c._fastqscreen_command() == (
                'path/to/fastqscreen --aligner bowtie2 r1.fastq.gz r2.fastq.gz '
                '--conf path/to/fastqscreen/conf --force'
            )

    def test_get_expected_outfiles(self):
        assert self.c.fastqscreen_expected_outfiles == 'r?_screen.txt'

    @patch('egcg_core.executor.execute')
    def test_run_fastqscreen(self, mocked_execute):
        with patch(ppath + 'ContaminationCheck._fastqscreen_command', return_value='a_cmd'):
            self.c._run()

        mocked_execute.assert_called_once_with(
            'a_cmd', working_dir='path/to/jobs/test_sample', mem=10, cpus=2, job_name='fastqscreen'
        )


class TestVerifyBamID(QCTester):
    def setUp(self):
        super().setUp()
        self.vbi = VerifyBamID(dataset=self.dataset, bam_file='test_bam_file.bam')

    @patch('egcg_core.executor.execute')
    def test_contamination_check(self, mocked_execute):
        with patch(ppath + 'util.find_file', new=self.fake_find_file):
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
