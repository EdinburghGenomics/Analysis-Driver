from .quality_control_base import QualityControl
from analysis_driver.config import default as cfg
from egcg_core import executor

class Relatedness(QualityControl):
    def __init__(self, dataset, working_dir, gVCF_files, reference, project_id):
        super().__init__(dataset, working_dir)
        self.gVCF_files = gVCF_files
        self.reference = reference
        self.project_id = project_id

    def get_gatk_genotype_gvcfs_command(self):
        gVCF_variants = " ". join(["--variant " + i for i in self.gVCF_files])
        out_prefix = self.project_id + '_genotype_gvcfs'
        number_threads = 12
        cmd = 'java -jar %s -T GenotypeGVCFs -nt %s -R %s %s -o %s' % (cfg['tools']['gatk'], number_threads, self.reference, gVCF_variants, out_prefix)
        gatk_outfile = out_prefix + '.vcf'
        return cmd, gatk_outfile

    def get_vcftools_relatedness_command(self, vcftools_input_file):
        cmd = '%s --relatedness2 --vcf %s' % (cfg['tools']['vcftools'], vcftools_input_file)
        vcftools_relatedness_outfile = vcftools_input_file.split('.')[0] + '.relatedness2'
        return cmd, vcftools_relatedness_outfile

    def run_gatk(self):
        gatk_genotype_gvcfs_command, gatk_outfile = self.get_gatk_genotype_gvcfs_command()
        self.dataset.start_stage('gatk_genotype_gvcfs')
        gatk_genotype_gvcfs_executor = executor.execute(
            gatk_genotype_gvcfs_command,
            job_name='gatk_genotype_gvcfs',
            working_dir=self.working_dir,
            cpus=12,
            mem=10
        )
        exit_status = gatk_genotype_gvcfs_executor.join()
        self.dataset.end_stage('gatk_genotype_gvcfs', exit_status)
        return gatk_outfile

    def run_vcftools(self, vcftools_input_file):
        vcftools_relatedness_command, vcftools_outfile = self.get_vcftools_relatedness_command(vcftools_input_file)
        self.dataset.start_stage('vcftools_relatedness')
        vcftools_relatedness_executor = executor.execute(
            vcftools_relatedness_command,
            job_name='vcftools_relatedness',
            working_dir=self.working_dir,
            cpus=1,
            mem=10
        )
        exit_status = vcftools_relatedness_executor.join()
        self.dataset.end_stage('vcftools_relatedness', exit_status)
        return vcftools_outfile

    def calculate_relatedness(self):
        gatk_genotype_gvcfs_outfile = self.run_gatk()
        vcftools_relatedness_outfile = self.run_vcftools(gatk_genotype_gvcfs_outfile)
        return vcftools_relatedness_outfile

    def run(self):
        try:
            self.vcftools_relatedness_expected_outfile = self.calculate_relatedness()
        except Exception as e:
            self.exception = e
            self.exit_status = 8

    def join(self, timeout=None):
        super().join(timeout=timeout)
        if self.exception:
            raise self.exception
        return self.vcftools_relatedness_expected_outfile