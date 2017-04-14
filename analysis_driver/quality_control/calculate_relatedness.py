from egcg_core import executor
from luigi import Parameter, ListParameter
from analysis_driver.config import default as cfg
from analysis_driver.segmentation import Stage


class Relatedness(Stage):
    gvcf_files = ListParameter()
    reference = Parameter()
    project_id = Parameter()

    @property
    def vcftools_relatedness_outfile(self):
        return self.project_id + '.relatedness2'

    @property
    def gatk_outfile(self):
        return self.project_id + '_genotype_gvcfs.vcf'

    def get_gatk_genotype_gvcfs_command(self):
        gvcf_variants = ' '. join(['--variant ' + i for i in self.gvcf_files])
        number_threads = 12
        return 'java -jar %s -T GenotypeGVCFs -nt %s -R %s %s -o %s' % (
            cfg['tools']['gatk'], number_threads, self.reference, gvcf_variants, self.gatk_outfile
        )

    def get_vcftools_relatedness_command(self):
        return '%s --relatedness2 --vcf %s --out %s' % (cfg['tools']['vcftools'], self.gatk_outfile, self.project_id)

    def run_gatk(self):
        gatk_genotype_gvcfs_command = self.get_gatk_genotype_gvcfs_command()
        return executor.execute(
            gatk_genotype_gvcfs_command,
            job_name='gatk_genotype_gvcfs',
            working_dir=self.job_dir,
            cpus=12,
            mem=30
        ).join()

    def run_vcftools(self):
        cmd = self.get_vcftools_relatedness_command()
        return executor.execute(
            cmd,
            job_name='vcftools_relatedness',
            working_dir=self.job_dir,
            cpus=1,
            mem=10
        ).join()

    def calculate_relatedness(self):
        exit_status = self.run_gatk()
        return exit_status + self.run_vcftools()
