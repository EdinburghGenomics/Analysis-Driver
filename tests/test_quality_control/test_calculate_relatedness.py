from unittest.mock import patch, Mock
from tests.test_quality_control.qc_tester import QCTester
from analysis_driver.quality_control.calculate_relatedness import Relatedness, Peddy, Genotype_gVCFs

ppath = 'analysis_driver.quality_control.calculate_relatedness.'


class TestGenotype_gVCFs(QCTester):
    def setUp(self):
        super().setUp()
        self.g = Genotype_gVCFs(
            dataset=self.project_dataset,
            gVCFs=['test_sample1.g.vcf.gz', 'test_sample2.g.vcf.gz', 'test_sample3.g.vcf.gz'],
            reference='/path/to/reference.fa'
        )

    def test_gatk_genotype_gvcfs_cmd(self):
        with patch(ppath + 'util.find_file', new=self.fake_find_file):
            assert self.r.gatk_genotype_gvcfs_cmd() == (
                'java -Djava.io.tmpdir=path/to/jobs/test_sample -XX:+UseSerialGC -Xmx20G '
                '-jar path/to/GenomeAnalysisTK.jar -T GenotypeGVCFs -nt 12 -R /path/to/reference.fa '
                '--variant test_sample1.g.vcf.gz --variant test_sample2.g.vcf.gz '
                '--variant test_sample3.g.vcf.gz -o test_project_id_genotype_gvcfs.vcf'
            )

    @patch('egcg_core.executor.execute')
    def test_run_gatk(self, mocked_execute):
        with patch(ppath + 'Genotype_gVCFs.gatk_genotype_gvcfs_cmd', return_value='test_command'):
            self.g.run_gatk()
            mocked_execute.assert_called_with(
                'test_command',
                job_name='gatk_genotype_gvcfs',
                cpus=12,
                working_dir='path/to/jobs/test_project_id',
                mem=30
            )


class TestRelatedness(QCTester):
    def setUp(self):
        super().setUp()
        self.r = Relatedness(
            dataset=self.project_dataset,
            genotyped_gvcfs='test_project_id_genotype_gvcfs.vcf',
            project_id='test_project_id'
        )

    def test_vcftools_relatedness_cmd(self):
        exp = 'path/to/vcftools --relatedness2 --vcf test_project_id_genotype_gvcfs.vcf --out test_project_id'
        assert self.r.vcftools_relatedness_cmd() == exp

    @patch('egcg_core.executor.execute')
    def test_run_vcftools(self, mocked_execute):
        with patch(ppath + 'Relatedness.vcftools_relatedness_cmd', return_value='test_command'):
            self.r.run_vcftools()
        mocked_execute.assert_called_with(
            'test_command',
            job_name='vcftools_relatedness',
            working_dir='path/to/jobs/test_project_id',
            cpus=1,
            mem=10
        )


class TestPeddy(QCTester):
    def setUp(self):
        super().setUp()
        self.p = Peddy(dataset=self.project_dataset,
                       genotyped_gvcfs='test_project_id_genotype_gvcfs.vcf',
                       ids=['test_sample1', 'test_sample2', 'test_sample3'],
                       )

    @patch('egcg_core.executor.execute')
    @patch(ppath + 'Peddy.peddy_command', return_value='peddy_command')
    def test_run_peddy(self, mocked_peddy_command, mocked_execute):
        self.p.run_peddy()
        mocked_execute.assert_called_with(
                mocked_peddy_command,
                job_name='peddy',
                cpus=10,
                working_dir='path/to/jobs/test_project_id',
                mem=10
            )

    @patch(ppath + 'Peddy.write_ped_file', return_value='ped.fam')
    def test_peddy_command(self, mocked_ped_file):
        assert self.p.peddy_command == 'peddy --plot --prefix test_project_id test_project_id_genotype_gvcfs.vcf ped.fam'

    @patch(ppath + 'Peddy.family_id', side_effect=['FAM1', 'FAM1', 'FAM1'])
    @patch(ppath + 'Peddy.relationship', side_effect=['Mother', 'Father', 'Proband'])
    @patch(ppath + 'Peddy.sex', side_effect=['Female', 'Male', 'Male'])
    def test_ped_file_content(self, mocked_sex, mocked_relationships, mocked_families):
        assert self.p.ped_file_content == [['FAM1', 'test_sample1', '0', '0', '2', '0'],
                                           ['FAM1', 'test_sample2', '0', '0', '1', '0'],
                                           ['FAM1', 'test_sample3', 'test_sample2', 'test_sample1', '1', '0']]

    @patch('egcg_core.executor.execute')
    def test_tabix_index(self, mocked_execute):
        self.p.tabix_index()
        mocked_execute.assert_called_with('tabix -f -p vcf test_project_id_genotype_gvcfs.vcf',
                                          job_name='tabix',
                                          cpus=12,
                                          mem=30,
                                          working_dir='path/to/jobs/test_project_id')

    def test_tabix_command(self):
        assert self.p.tabix_command == 'tabix -f -p vcf test_project_id_genotype_gvcfs.vcf'
