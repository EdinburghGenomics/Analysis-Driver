from unittest.mock import patch
from tests.test_quality_control.qc_tester import QCTester
from analysis_driver.quality_control.relatedness import Relatedness

ppath = 'analysis_driver.quality_control.relatedness.'


class TestRelatedness(QCTester):
    def setUp(self):
        super().setUp()
        self.r = Relatedness(
            dataset=self.dataset,
            gvcf_files=['test_sample1.g.vcf.gz', 'test_sample2.g.vcf.gz', 'test_sample3.g.vcf.gz'],
            reference='/path/to/reference.fa',
            project_id='test_project_id'
        )

    def test_gatk_genotype_gvcfs_cmd(self):
        with patch(ppath + 'util.find_file', new=self.fake_find_file):
            assert self.r.gatk_genotype_gvcfs_cmd() == (
                'java -Djava.io.tmpdir=path/to/jobs/test_sample -XX:+UseSerialGC -Xmx20G '
                '-jar path/to/GenomeAnalysisTK.jar -T GenotypeGVCFs -nt 12 -R /path/to/reference.fa '
                '--variant test_sample1.g.vcf.gz --variant test_sample2.g.vcf.gz '
                '--variant test_sample3.g.vcf.gz -o test_project_id_genotype_gvcfs.vcf'
            )

    def test_vcftools_relatedness_cmd(self):
        exp = 'path/to/vcftools --relatedness2 --vcf test_project_id_genotype_gvcfs.vcf --out test_project_id'
        assert self.r.vcftools_relatedness_cmd() == exp

    @patch('egcg_core.executor.execute')
    def test_run_gatk(self, mocked_execute):
        with patch(ppath + 'Relatedness.gatk_genotype_gvcfs_cmd', return_value='test_command'):
            self.r.run_gatk()
            mocked_execute.assert_called_with(
                'test_command',
                job_name='gatk_genotype_gvcfs',
                cpus=12,
                working_dir='path/to/jobs/test_sample',
                mem=30
            )

    @patch('egcg_core.executor.execute')
    def test_run_vcftools(self, mocked_execute):
        with patch(ppath + 'Relatedness.vcftools_relatedness_cmd', return_value='test_command'):
            self.r.run_vcftools()
        mocked_execute.assert_called_with(
            'test_command',
            job_name='vcftools_relatedness',
            working_dir='path/to/jobs/test_sample',
            cpus=1,
            mem=10
        )
