from unittest.mock import patch
from tests.test_quality_control.qc_tester import QCTester
from analysis_driver.quality_control.calculate_relatedness import Relatedness

class Test_Relatedness(QCTester):
    def setUp(self):
        super().setUp()
        self.working_dir = 'test_project'
        self.gVCF_files = ['test_sample1.g.vcf.gz',
                           'test_sample2.g.vcf.gz',
                           'test_sample3.g.vcf.gz']
        self.reference = '/path/to/reference.fa'
        self.project_id = 'test_project_id'
        self.r = Relatedness(self.dataset, self.working_dir, self.gVCF_files, self.reference, self.project_id)

    def test_get_gatk_genotype_gvcfs_command(self):
        cmd, outfile = self.r.get_gatk_genotype_gvcfs_command()
        assert cmd == 'java -jar path/to/GenomeAnalysisTK.jar ' \
                      '-T GenotypeGVCFs ' \
                      '-R /path/to/reference.fa ' \
                      '--variant test_sample1.g.vcf.gz ' \
                      '--variant test_sample2.g.vcf.gz ' \
                      '--variant test_sample3.g.vcf.gz ' \
                      '-o test_project_id_genotype_gvcfs'
        assert outfile == 'test_project_id_genotype_gvcfs.vcf'

    def test_get_vcftools_relatedness_command(self):
        input_vcf = 'test_project_id_genotype_gvcfs.vcf'
        cmd, vcftools_outfile = self.r.get_vcftools_relatedness_command(input_vcf)
        assert cmd == 'path/to/vcftools --relatedness --vcf test_project_id_genotype_gvcfs.vcf'
        assert vcftools_outfile == 'test_project_id_genotype_gvcfs.relatedness2'

    @patch('analysis_driver.quality_control.calculate_relatedness.Relatedness.get_gatk_genotype_gvcfs_command')
    @patch('analysis_driver.dataset.rest_communication')
    @patch('egcg_core.executor.execute')
    def test_run_gatk(self, mocked_execute, mocked_rest, mocked_gatk_outfile):
        mocked_gatk_outfile.return_value = ('test_command', 'test_outfile')
        run_gatk_outfile = self.r.run_gatk()
        mocked_execute.assert_called_once_with('test_command', job_name='gatk_genotype_gvcfs', cpus=2, working_dir='test_project', mem=10)
        assert run_gatk_outfile == 'test_outfile'

    @patch('analysis_driver.quality_control.calculate_relatedness.Relatedness.get_vcftools_relatedness_command')
    @patch('analysis_driver.dataset.rest_communication')
    @patch('egcg_core.executor.execute')
    def test_run_vcftools(self, mocked_execute, mocked_rest, mocked_vcftools_outfile):
        mocked_vcftools_outfile.return_value = ('test_command', 'test_outfile')
        run_vcftools_outfile = self.r.run_vcftools('test_input_file')
        mocked_execute.assert_called_once_with('test_command', job_name='vcftools_relatedness', working_dir='test_project', cpus=2, mem=10)