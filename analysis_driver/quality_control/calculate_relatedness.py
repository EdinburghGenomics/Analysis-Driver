from egcg_core import executor, util
from luigi import Parameter, ListParameter
from analysis_driver.config import default as cfg
from analysis_driver.segmentation import Stage


class Relatedness(Stage):
    gvcf_files = ListParameter()
    reference = Parameter()
    project_id = Parameter()

    @property
    def gatk_outfile(self):
        return self.project_id + '_genotype_gvcfs.vcf'

    def gatk_genotype_gvcfs_cmd(self):
        gvcf_variants = ' '. join(['--variant ' + util.find_file(i) for i in self.gvcf_files])
        number_threads = 12
        return 'java -jar %s -T GenotypeGVCFs -nt %s -R %s %s -o %s' % (
            cfg['tools']['gatk'], number_threads, self.reference, gvcf_variants, self.gatk_outfile
        )

    def vcftools_relatedness_cmd(self):
        return '%s --relatedness2 --vcf %s --out %s' % (
            cfg['tools']['vcftools'], self.gatk_outfile, self.project_id
        )

    def run_gatk(self):
        return executor.execute(
            self.gatk_genotype_gvcfs_cmd(),
            job_name='gatk_genotype_gvcfs',
            working_dir=self.job_dir,
            cpus=12,
            mem=30
        ).join()

    def run_vcftools(self):
        return executor.execute(
            self.vcftools_relatedness_cmd(),
            job_name='vcftools_relatedness',
            working_dir=self.job_dir,
            cpus=1,
            mem=10
        ).join()

    def _run(self):
        return self.run_gatk() + self.run_vcftools()
