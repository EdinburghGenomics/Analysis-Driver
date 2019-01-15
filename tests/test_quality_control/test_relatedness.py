import os
from unittest.mock import patch, Mock
from tests.test_quality_control.qc_tester import QCTester
from analysis_driver.quality_control.relatedness import Relatedness, GenotypeGVCFs, Peddy, ParseRelatedness
from analysis_driver.exceptions import PipelineError


ppath = 'analysis_driver.quality_control.relatedness.'


class TestGenotypeGVCFs(QCTester):
    def setUp(self):
        super().setUp()
        self.g = GenotypeGVCFs(
            dataset=self.project_dataset,
            gVCFs=['test_sample1.g.vcf.gz', 'test_sample2.g.vcf.gz', 'test_sample3.g.vcf.gz'],
            reference='/path/to/reference.fa'
        )

    def test_memory(self):
        assert self.g.memory == 50

    def test_gatk_genotype_gvcfs_cmd(self):
        with patch(ppath + 'util.find_file', new=self.fake_find_file):
            assert self.g.gatk_genotype_gvcfs_cmd(50, number_threads=12) == (
                'java -Djava.io.tmpdir=tests/assets/jobs/test_project_id -XX:+UseSerialGC -Xmx50G '
                '-jar path/to/gatk -T GenotypeGVCFs -nt 12 -R /path/to/reference.fa '
                '--variant test_sample1.g.vcf.gz --variant test_sample2.g.vcf.gz '
                '--variant test_sample3.g.vcf.gz -o tests/assets/jobs/test_project_id/test_project_id_genotype_gvcfs.vcf'
            )

    @patch('egcg_core.executor.execute')
    def test_run(self, mocked_execute):
        with patch(ppath + 'GenotypeGVCFs.gatk_genotype_gvcfs_cmd', return_value='test_command'):
            self.g._run()
            mocked_execute.assert_called_with(
                'test_command',
                job_name='gatk_genotype_gvcfs',
                cpus=12,
                working_dir='tests/assets/jobs/test_project_id',
                mem=52
            )


class TestRelatedness(QCTester):
    def setUp(self):
        super().setUp()
        self.r = Relatedness(dataset=self.project_dataset)

    def test_vcftools_relatedness_cmd(self):
        exp = 'path/to/vcftools --relatedness2 --vcf tests/assets/jobs/test_project_id/test_project_id_genotype_gvcfs.vcf --out test_project_id'
        assert self.r.vcftools_relatedness_cmd() == exp

    @patch('egcg_core.executor.execute')
    def test_run(self, mocked_execute):
        with patch(ppath + 'Relatedness.vcftools_relatedness_cmd', return_value='test_command'):
            self.r._run()
        mocked_execute.assert_called_with(
            'test_command',
            job_name='vcftools_relatedness',
            working_dir='tests/assets/jobs/test_project_id',
            cpus=1,
            mem=10
        )


class TestPeddy(QCTester):
    def setUp(self):
        super().setUp()
        self.p = Peddy(dataset=self.project_dataset, ids=['test_sample1', 'test_sample2', 'test_sample3'])

    @patch('egcg_core.executor.execute')
    @patch(ppath + 'Peddy.peddy_command', return_value='peddy_command')
    @patch(ppath + 'Peddy.write_ped_file')
    def test_run_peddy(self, mocked_write_ped, mocked_peddy_command, mocked_execute):
        self.p.run_peddy()
        mocked_execute.assert_called_with(
            mocked_peddy_command,
            job_name='peddy',
            cpus=1,
            working_dir='tests/assets/jobs/test_project_id',
            mem=10
        )

    def test_peddy_command(self):
        assert self.p.peddy_command == ('path/to/peddy --plot --prefix test_project_id '
                                        'tests/assets/jobs/test_project_id/test_project_id_genotype_gvcfs.vcf.gz '
                                        'tests/assets/jobs/test_project_id/ped.fam')

    @patch(ppath + 'Peddy.relationship')
    def test_relationships(self, mocked_relationship):
        family = ['test_sample1', 'test_sample2', 'test_sample3']

        mocked_relationship.side_effect = ['Father', 'Mother', 'Proband']
        assert self.p.relationships(family) == {'Proband': {'Mother': 'test_sample2', 'Father': 'test_sample1'},
                                                'Mother': {'Mother': '0', 'Father': '0'},
                                                'Father': {'Mother': '0', 'Father': '0'},
                                                'Sister': {'Mother': 'test_sample2', 'Father': 'test_sample1'},
                                                'Brother': {'Mother': 'test_sample2', 'Father': 'test_sample1'},
                                                'Other': {'Mother': '0', 'Father': '0'}}

        mocked_relationship.side_effect = ['Proband', 'Father', 'Mother']
        assert self.p.relationships(family) == {'Proband': {'Mother': 'test_sample3', 'Father': 'test_sample2'},
                                                'Mother': {'Mother': '0', 'Father': '0'},
                                                'Father': {'Mother': '0', 'Father': '0'},
                                                'Sister': {'Mother': 'test_sample3', 'Father': 'test_sample2'},
                                                'Brother': {'Mother': 'test_sample3', 'Father': 'test_sample2'},
                                                'Other': {'Mother': '0', 'Father': '0'}}

        mocked_relationship.side_effect = ['Proband', 'Other', 'Other']
        assert self.p.relationships(family) == {'Proband': {'Mother': '0', 'Father': '0'},
                                                'Mother': {'Mother': '0', 'Father': '0'},
                                                'Father': {'Mother': '0', 'Father': '0'},
                                                'Sister': {'Mother': '0', 'Father': '0'},
                                                'Brother': {'Mother': '0', 'Father': '0'},
                                                'Other': {'Mother': '0', 'Father': '0'}}

        family = ['test_sample1', 'test_sample2']
        mocked_relationship.side_effect = ['Mother', 'Proband']
        assert self.p.relationships(family) == {'Proband': {'Mother': 'test_sample1', 'Father': '0'},
                                                'Mother': {'Mother': '0', 'Father': '0'},
                                                'Father': {'Mother': '0', 'Father': '0'},
                                                'Sister': {'Mother': 'test_sample1', 'Father': '0'},
                                                'Brother': {'Mother': 'test_sample1', 'Father': '0'},
                                                'Other': {'Mother': '0', 'Father': '0'}}

    @patch(ppath + 'Peddy.family_id')
    @patch(ppath + 'Peddy.get_member_details')
    @patch(ppath + 'Peddy.all_families')
    def test_ped_file_content(self, pfams, pmem, pfam):
        pfam.side_effect = ['FAM1', 'FAM1', 'FAM2']
        pmem.side_effect = [[['FAM1', '0', '0', '2', '0'], ['FAM1', 'test_sample1', '0', '1', '0']], [['FAM2', '0', '0', '1', '0']]]
        pfams.return_value = {'FAM1': ['test_sample1', 'test_sample2'], 'FAM2': ['test_sample3']}
        assert self.p.ped_file_content == [['FAM1', '0', '0', '2', '0'], ['FAM1', 'test_sample1', '0', '1', '0'], ['FAM2', '0', '0', '1', '0']]

        pfam.side_effect = ['FAM1', 'FAM1', 'FAM1']
        pmem.side_effect = [[['FAM1', '0', '0', '2', '0'], ['FAM1', '0', '0', '1', '0'], ['FAM1', 'test_sample1', 'test_sample2', '1', '0']]]
        pfams.return_value = {'FAM1': ['test_sample1', 'test_sample2', 'test_sample3']}
        assert self.p.ped_file_content == [['FAM1', '0', '0', '2', '0'], ['FAM1', '0', '0', '1', '0'], ['FAM1', 'test_sample1', 'test_sample2', '1', '0']]

        pfam.side_effect = ['FAM1', 'FAM1', 'No_ID']
        pmem.side_effect = [[['FAM1', '0', '0', '2', '0'], ['FAM1', '0', 'test_sample1', '1', '0']], [['No_ID', '0', '0', '1', '0']]]
        pfams.return_value = {'FAM1': ['test_sample1', 'test_sample2'], 'No_ID': ['test_sample3']}
        assert self.p.ped_file_content == [['FAM1', '0', '0', '2', '0'], ['FAM1', '0', 'test_sample1', '1', '0'], ['No_ID', '0', '0', '1', '0']]

    @patch(ppath + 'Peddy.relationships')
    @patch(ppath + 'Peddy.relationship')
    @patch(ppath + 'sex_check.alias')
    @patch(ppath + 'clarity.get_user_sample_name')
    @patch(ppath + 'clarity.get_sample')
    def test_get_member_details(self, psample, pname, psex, prel, prels):
        psample.return_value = Mock(udf={'Sex': 'F'})
        prels.return_value = {'Proband': {'Mother': 'test_sample1', 'Father': '0'},
                              'Mother': {'Mother': '0', 'Father': '0'},
                              'Father': {'Mother': '0', 'Father': '0'},
                              'Sister': {'Mother': 'test_sample1', 'Father': '0'},
                              'Brother': {'Mother': 'test_sample1', 'Father': '0'},
                              'Other': {'Mother': '0', 'Father': '0'}}
        prel.side_effect = ['Mother', 'Proband']
        psex.side_effect = ['female', 'male']
        pname.side_effect = ['usersample1', 'usersample2', 'usersample1']
        all_families = {'FAM1': ['test_sample1', 'test_sample2'], 'FAM2': ['test_sample3']}
        assert self.p.get_member_details('FAM1', all_families) == [
            ['FAM1', 'usersample1', '0', '0', '2', '0'], ['FAM1', 'usersample2', '0', 'usersample1', '1', '0']
        ]

        prels.return_value = {'Proband': {'Mother': '0', 'Father': '0'},
                              'Mother': {'Mother': '0', 'Father': '0'},
                              'Father': {'Mother': '0', 'Father': '0'},
                              'Sister': {'Mother': '0', 'Father': '0'},
                              'Brother': {'Mother': '0', 'Father': '0'},
                              'Other': {'Mother': '0', 'Father': '0'}}

        prel.side_effect = ['Other', 'Other', 'Proband']
        psex.side_effect = ['unknown', 'unknown', 'male']
        pname.side_effect = ['usersample1', 'usersample2', 'usersample3']
        all_families = {'FAM1': ['test_sample1', 'test_sample2', 'test_sample3']}
        assert self.p.get_member_details('FAM1', all_families) == [['FAM1', 'usersample1', '0', '0', '0', '0'],
                                                                   ['FAM1', 'usersample2', '0', '0', '0', '0'],
                                                                   ['FAM1', 'usersample3', '0', '0', '1', '0']]

    @patch('egcg_core.executor.execute')
    def test_tabix_index(self, mocked_execute):
        self.p.tabix_index()
        mocked_execute.assert_called_with(
            'path/to/tabix -f -p vcf tests/assets/jobs/test_project_id/test_project_id_genotype_gvcfs.vcf.gz',
            job_name='tabix',
            cpus=1,
            mem=10,
            working_dir='tests/assets/jobs/test_project_id'
        )

    @patch('egcg_core.executor.execute')
    def test_bgzip(self, mocked_execute):
        self.p.bgzip()
        mocked_execute.assert_called_with(
            'path/to/bgzip -f tests/assets/jobs/test_project_id/test_project_id_genotype_gvcfs.vcf',
            job_name='bgzip',
            cpus=1,
            mem=10,
            working_dir='tests/assets/jobs/test_project_id'
        )


class TestParseRelatedness(QCTester):
    def setUp(self):
        super().setUp()
        self.p = ParseRelatedness(dataset=self.project_dataset, parse_method='parse_both', ids=['test_sample1', 'test_sample2', 'test_sample3'])

    @patch(ppath + 'clarity.get_user_sample_name')
    def test_user_sample_ids(self, pname):
        pname.side_effect = ['user_sample1', 'user_sample2', 'user_sample3']
        assert self.p.user_sample_ids() == {
            'user_sample1': 'test_sample1', 'user_sample2': 'test_sample2', 'user_sample3': 'test_sample3'
        }

        pname.side_effect = ['user_sample1', 'user_sample1', 'user_sample2']
        with self.assertRaises(PipelineError) as c:
            self.p.user_sample_ids()
        assert str(c.exception) == 'User ID user_sample1 appears more than once in sample list'

    @patch(ppath + 'ParseRelatedness.user_sample_ids')
    @patch(ppath + 'ParseRelatedness.family_id')
    @patch(ppath + 'ParseRelatedness.relationship')
    def test_get_outfile_content(self, prel, pfam, pname):
        pname.return_value = {
            'user_sample1': 'test_sample1',
            'user_sample2': 'test_sample2',
            'user_sample3': 'test_sample3',
            'user_sample4': 'test_sample4'
        }
        pfam.side_effect = ['FAM1', 'FAM1', 'FAM1', 'FAM1', 'FAM1', 'FAM1', 'FAM1', 'FAM1']
        prel.side_effect = ['Other', 'Other', 'Other', 'Other', 'Other', 'Other', 'Other', 'Other']

        gel, egc = self.p.get_outfile_content(
            [
                {'sample1': 'user_sample1', 'sample2': 'user_sample2', 'relatedness': [1, 0.9]},
                {'sample1': 'user_sample3', 'sample2': 'user_sample4', 'relatedness': [0.9, 0.7]}
            ]
        )
        assert gel == []
        assert egc == [['test_sample1', 'FAM1', 'Other', 'test_sample2', 'FAM1', 'Other', 1, 0.9],
                       ['test_sample3', 'FAM1', 'Other', 'test_sample4', 'FAM1', 'Other', 0.9, 0.7]]

        pfam.side_effect = ['FAM1', 'FAM1', 'FAM1', 'FAM1']
        prel.side_effect = ['Proband', 'Other', 'Other', 'Other']
        pname.return_value = {'user_sample1': 'test_sample1', 'user_sample2': 'test_sample2',
                              'user_sample3': 'test_sample3', 'user_sample4': 'test_sample4'}

        gel, egc = self.p.get_outfile_content(
            [{'sample1': 'user_sample1', 'sample2': 'user_sample2', 'relatedness': [1, 0.9]},
             {'sample1': 'user_sample3', 'sample2': 'user_sample4', 'relatedness': [0.9, 0.7]}]
        )
        assert gel == [['FAM1', 'user_sample1', 'Proband', 'user_sample2', 'Other', 0.9]]
        assert egc == [['test_sample1', 'FAM1', 'Proband', 'test_sample2', 'FAM1', 'Other', 1, 0.9],
                       ['test_sample3', 'FAM1', 'Other', 'test_sample4', 'FAM1', 'Other', 0.9, 0.7]]

        pfam.side_effect = ['FAM1', 'FAM1', 'FAM1', 'FAM1', 'FAM1', 'FAM1', 'FAM1', 'FAM1']
        prel.side_effect = ['Proband', 'Other', 'Other', 'Proband', 'Other', 'Other', 'Other', 'Other', 'Other']
        pname.return_value = {'user_sample1': 'test_sample1', 'user_sample2': 'test_sample2',
                              'user_sample3': 'test_sample3', 'user_sample4': 'test_sample4'}

        gel, egc = self.p.get_outfile_content(
            [{'sample1': 'user_sample1', 'sample2': 'user_sample1', 'relatedness': [1, 0.9]},
             {'sample1': 'user_sample2', 'sample2': 'user_sample2', 'relatedness': [0.9, 0.7]}]
        )
        assert gel == []
        assert egc == []

    def test_get_columns(self):
        peddy_file = os.path.join(self.assets_path, self.project_id + '.ped_check.csv')
        assert self.p.get_columns(peddy_file, ['sample_a', 'sample_b', 'rel']) == [
            ['NA12777', 'NA12777NA12877', '1'],
            ['NA12777', 'NA12877', '1'],
            ['NA12777', 'NA12878', '0']
        ]

        with self.assertRaises(FileNotFoundError):
            self.p.get_columns('/path/to/nonexistent/file', ['sample_a', 'sample_b', 'rel'])

    @patch(ppath + 'ParseRelatedness.relatedness_file', new=os.path.join(QCTester.assets_path, 'test_project.relatedness2'))
    @patch(ppath + 'ParseRelatedness.get_columns')
    def test_combine_peddy_vcftools(self, mocked_get_cols):
        mocked_get_cols.return_value = [
            ['NA12777', 'NA12777NA12877', '1'],
            ['NA12777', 'NA12877', '-0.09722'],
            ['NA12777', 'NA12878', '-0.07692'],
            ['NA12777', 'NA12882', '-0.1667'],
            ['NA12777NA12877', 'NA12877', '0.8077'],
            ['NA12777NA12877', 'NA12878', '0.4923'],
            ['NA12777NA12877', 'NA12882', '0.6364'],
            ['NA12877', 'NA12878', '-0.03077'],
            ['NA12877', 'NA12882', '0.5303'],
            ['NA12878', 'NA12882', '0.4769']
        ]
        assert self.p.combine_peddy_vcftools() == [
            {'sample2': 'NA12777NA12877', 'sample1': 'NA12777', 'relatedness': ['1', '0.369531']},
            {'sample2': 'NA12877', 'sample1': 'NA12777', 'relatedness': ['-0.09722', '-0.00872645']},
            {'sample2': 'NA12878', 'sample1': 'NA12777', 'relatedness': ['-0.07692', '-0.00330048']},
            {'sample2': 'NA12882', 'sample1': 'NA12777', 'relatedness': ['-0.1667', '-0.00613699']},
            {'sample2': 'NA12877', 'sample1': 'NA12777NA12877', 'relatedness': ['0.8077', '0.28758']},
            {'sample2': 'NA12878', 'sample1': 'NA12777NA12877', 'relatedness': ['0.4923', '0.150708']},
            {'sample2': 'NA12882', 'sample1': 'NA12777NA12877', 'relatedness': ['0.6364', '0.213526']},
            {'sample2': 'NA12878', 'sample1': 'NA12877', 'relatedness': ['-0.03077', '0.00116826']},
            {'sample2': 'NA12882', 'sample1': 'NA12877', 'relatedness': ['0.5303', '0.218896']},
            {'sample2': 'NA12882', 'sample1': 'NA12878', 'relatedness': ['0.4769', '0.243608']}
        ]
