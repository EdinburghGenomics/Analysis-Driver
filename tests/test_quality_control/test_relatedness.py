from unittest.mock import patch
from tests.test_quality_control.qc_tester import QCTester
from analysis_driver.quality_control.relatedness import Relatedness, Genotype_gVCFs, Peddy

ppath = 'analysis_driver.quality_control.relatedness.'


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
            assert self.g.gatk_genotype_gvcfs_cmd() == (
                'java -Djava.io.tmpdir=path/to/jobs/test_project_id -XX:+UseSerialGC -Xmx50G '
                '-jar path/to/gatk -T GenotypeGVCFs -nt 12 -R /path/to/reference.fa '
                '--variant test_sample1.g.vcf.gz --variant test_sample2.g.vcf.gz '
                '--variant test_sample3.g.vcf.gz -o path/to/jobs/test_project_id/test_project_id_genotype_gvcfs.vcf'
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
                mem=50
            )


class TestRelatedness(QCTester):
    def setUp(self):
        super().setUp()
        self.r = Relatedness(
            dataset=self.project_dataset,
        )

    def test_vcftools_relatedness_cmd(self):
        exp = 'path/to/vcftools --relatedness2 --vcf path/to/jobs/test_project_id/test_project_id_genotype_gvcfs.vcf --out test_project_id'
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
                       ids=['test_sample1', 'test_sample2', 'test_sample3'],
                       )

    @patch('egcg_core.executor.execute')
    @patch(ppath + 'Peddy.peddy_command', return_value='peddy_command')
    def test_run_peddy(self, mocked_peddy_command, mocked_execute):
        self.p.run_peddy()
        mocked_execute.assert_called_with(
                mocked_peddy_command,
                job_name='peddy',
                cpus=1,
                working_dir='path/to/jobs/test_project_id',
                mem=10
            )

    @patch(ppath + 'Peddy.write_ped_file', return_value='ped.fam')
    def test_peddy_command(self, mocked_ped_file):
        assert self.p.peddy_command == 'path/to/peddy --plot --prefix test_project_id path/to/jobs/test_project_id/test_project_id_genotype_gvcfs.vcf.gz ped.fam'

    def test_relationships(self):
        family = ['test_sample1', 'test_sample2', 'test_sample3']

        with patch(ppath + 'Peddy.relationship', side_effect=['Father', 'Mother', 'Proband']):
            assert self.p.relationships(family) == {'Proband': {'Mother':'test_sample2', 'Father':'test_sample1'},
                                                      'Mother':{'Mother':'0', 'Father':'0'},
                                                      'Father':{'Mother':'0', 'Father':'0'},
                                                      'Sister':{'Mother':'test_sample2', 'Father':'test_sample1'},
                                                      'Brother': {'Mother':'test_sample2', 'Father':'test_sample1'},
                                                      'Other':{'Mother':'0', 'Father':'0'}}

        with patch(ppath + 'Peddy.relationship', side_effect=['Proband', 'Father', 'Mother']):
            assert self.p.relationships(family) == {'Proband': {'Mother':'test_sample3', 'Father':'test_sample2'},
                                                      'Mother':{'Mother':'0', 'Father':'0'},
                                                      'Father':{'Mother':'0', 'Father':'0'},
                                                      'Sister':{'Mother':'test_sample3', 'Father':'test_sample2'},
                                                      'Brother': {'Mother':'test_sample3', 'Father':'test_sample2'},
                                                      'Other':{'Mother':'0', 'Father':'0'}}

        with patch(ppath + 'Peddy.relationship', side_effect=['Proband', 'Other', 'Other']):
            assert self.p.relationships(family) == {'Proband': {'Mother':'0', 'Father':'0'},
                                                      'Mother':{'Mother':'0', 'Father':'0'},
                                                      'Father':{'Mother':'0', 'Father':'0'},
                                                      'Sister':{'Mother':'0', 'Father':'0'},
                                                      'Brother': {'Mother':'0', 'Father':'0'},
                                                      'Other':{'Mother':'0', 'Father':'0'}}

        family = ['test_sample1', 'test_sample2']
        with patch(ppath + 'Peddy.relationship', side_effect=['Mother', 'Proband']):
            assert self.p.relationships(family) == {'Proband': {'Mother':'test_sample1', 'Father':'0'},
                                                      'Mother':{'Mother':'0', 'Father':'0'},
                                                      'Father':{'Mother':'0', 'Father':'0'},
                                                      'Sister':{'Mother':'test_sample1', 'Father':'0'},
                                                      'Brother': {'Mother':'test_sample1', 'Father':'0'},
                                                      'Other':{'Mother':'0', 'Father':'0'}}

    def test_ped_file_content(self):
        with patch(ppath + 'Peddy.family_id', side_effect=['FAM1', 'FAM1', 'FAM2']):
            with patch(ppath + 'Peddy.get_member_details', side_effect=[[['FAM1', '0', '0', '2', '0'], ['FAM1', 'test_sample1', '0', '1', '0']], [['FAM2', '0', '0', '1', '0']]]):
                with patch(ppath + 'Peddy.all_families', return_value={'FAM1': ['test_sample1', 'test_sample2'], 'FAM2': ['test_sample3']}):
                    assert self.p.ped_file_content == [['FAM1', '0', '0', '2', '0'], ['FAM1', 'test_sample1', '0', '1', '0'], ['FAM2', '0', '0', '1', '0']]

        with patch(ppath + 'Peddy.family_id', side_effect=['FAM1', 'FAM1', 'FAM1']):
            with patch(ppath + 'Peddy.get_member_details', side_effect=[[['FAM1', '0', '0', '2', '0'], ['FAM1', '0', '0', '1', '0'], ['FAM1', 'test_sample1', 'test_sample2', '1', '0']]]):
                with patch(ppath + 'Peddy.all_families', return_value={'FAM1': ['test_sample1', 'test_sample2', 'test_sample3']}):
                    assert self.p.ped_file_content == [['FAM1', '0', '0', '2', '0'], ['FAM1', '0', '0', '1', '0'], ['FAM1', 'test_sample1', 'test_sample2', '1', '0']]

        with patch(ppath + 'Peddy.family_id', side_effect=['FAM1', 'FAM1', 'No_ID']):
            with patch(ppath + 'Peddy.get_member_details', side_effect=[[['FAM1', '0', '0', '2', '0'], ['FAM1', '0', 'test_sample1', '1', '0']], [['No_ID', '0', '0', '1', '0']]]):
                with patch(ppath + 'Peddy.all_families', return_value={'FAM1': ['test_sample1', 'test_sample2'], 'No_ID': ['test_sample3']}):
                    assert self.p.ped_file_content == [['FAM1', '0', '0', '2', '0'],
                                                       ['FAM1', '0', 'test_sample1', '1', '0'],
                                                       ['No_ID', '0', '0', '1', '0']]

    def test_get_member_details(self):
        with patch(ppath + 'Peddy.relationships', return_value={'Proband': {'Mother':'test_sample1', 'Father':'0'},
                                                          'Mother':{'Mother':'0', 'Father':'0'},
                                                          'Father':{'Mother':'0', 'Father':'0'},
                                                          'Sister':{'Mother':'test_sample1', 'Father':'0'},
                                                          'Brother': {'Mother':'test_sample1', 'Father':'0'},
                                                          'Other':{'Mother':'0', 'Father':'0'}}):
            with patch(ppath + 'Peddy.relationship', side_effect=['Mother', 'Proband']):
                with patch(ppath + 'Peddy.sex', side_effect=['Female', 'Male']):
                    with patch(ppath + 'clarity.get_user_sample_name', side_effect=['usersample1', 'usersample2', 'usersample1']):
                        all_families = {'FAM1': ['test_sample1', 'test_sample2'], 'FAM2': ['test_sample3']}
                        assert self.p.get_member_details('FAM1', all_families) == [['FAM1', 'usersample1', '0', '0', '2', '0'],
                                                                                   ['FAM1', 'usersample2', '0', 'usersample1', '1', '0']]



        with patch(ppath + 'Peddy.relationships', return_value={'Proband': {'Mother':'0', 'Father':'0'},
                                                          'Mother':{'Mother':'0', 'Father':'0'},
                                                          'Father':{'Mother':'0', 'Father':'0'},
                                                          'Sister':{'Mother':'0', 'Father':'0'},
                                                          'Brother': {'Mother':'0', 'Father':'0'},
                                                          'Other':{'Mother':'0', 'Father':'0'}}):

            with patch(ppath + 'Peddy.relationship', side_effect=['Other', 'Other', 'Proband']):
                with patch(ppath + 'Peddy.sex', side_effect=['No_Sex', 'No_Sex', 'Male']):
                    with patch(ppath + 'clarity.get_user_sample_name', side_effect=['usersample1', 'usersample2', 'usersample3']):
                        all_families = {'FAM1': ['test_sample1', 'test_sample2', 'test_sample3']}
                        assert self.p.get_member_details('FAM1', all_families) == [['FAM1', 'usersample1', '0', '0', '0', '0'],
                                                                                   ['FAM1', 'usersample2', '0', '0', '0', '0'],
                                                                                   ['FAM1', 'usersample3', '0', '0', '1', '0']]

    @patch('egcg_core.executor.execute')
    def test_tabix_index(self, mocked_execute):
        self.p.tabix_index()
        mocked_execute.assert_called_with('path/to/tabix -f -p vcf path/to/jobs/test_project_id/test_project_id_genotype_gvcfs.vcf.gz',
                                          job_name='tabix',
                                          cpus=1,
                                          mem=10,
                                          working_dir='path/to/jobs/test_project_id')

    def test_tabix_command(self):
        assert self.p.tabix_command == 'path/to/tabix -f -p vcf path/to/jobs/test_project_id/test_project_id_genotype_gvcfs.vcf.gz'
