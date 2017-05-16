from egcg_core import executor, util
from luigi import Parameter, ListParameter
from analysis_driver.config import default as cfg
from analysis_driver.segmentation import Stage
from egcg_core import clarity

class Relatedness(Stage):
    gvcf_files = ListParameter()
    reference = Parameter()
    project_id = Parameter()

    @property
    def gatk_outfile(self):
        return self.dataset.name + '_genotype_gvcfs.vcf'

    def gatk_genotype_gvcfs_cmd(self, gvcf_files):
        gvcf_variants = ' '. join(['--variant ' + util.find_file(i) for i in gvcf_files])
        number_threads = 12
        return 'java -jar %s -T GenotypeGVCFs -nt %s -R %s %s -o %s' % (
            cfg['tools']['gatk'], number_threads, self.reference, gvcf_variants, self.gatk_outfile
        )

    def vcftools_relatedness_cmd(self):
        return '%s --relatedness2 --vcf %s --out %s' % (
            cfg['tools']['vcftools'], self.gatk_outfile, self.project_id
        )

    def run_gatk(self, gvcf_files):
        return executor.execute(
            self.gatk_genotype_gvcfs_cmd(gvcf_files),
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
        return self.run_gatk(self.gvcf_files) + self.run_vcftools()


class Peddy(Relatedness):
    gvcf_files = ListParameter()
    reference = Parameter()

    @property
    def tabix_command(self):
         return "tabix -f -p vcf %s" % (self.gatk_outfile)


    def tabix_index(self):
        return executor.execute(
            self.tabix_command,
            job_name='tabix',
            working_dir=self.job_dir,
            cpus=12,
            mem=30
        ).join()

    @property
    def ped_file(self):
        ped_file = 'ped.fam'
        sex_codes = {'Male': 1, 'Female': 2}
        all_families = {}
        for i in self.ids:
            family_id = clarity.get_sample(id).udf.get('Family ID')

            if not all_families[family_id]:
                all_families[family_id] = []
            all_families[family_id].append(i)

        with open(ped_file, 'w') as openfile:
            for family in all_families:
                relationship_codes = {'Proband': {'mother':None, 'father':None},
                                      'Mother':{'mother':0, 'father':0},
                                      'Father':{'mother':0, 'father':0},
                                      'Sister':{'mother':None, 'father':None},
                                      'Brother': {'mother':None, 'father':None},
                                      'Other':{'mother':0, 'father':0}}

                for member in all_families[family]:
                    relationship = clarity.get_sample(member).udf.get('Relationship')
                    if relationship == 'Father':
                        relationship_codes['Proband']['father'] = (member)
                        relationship_codes['Sister']['father'] = (member)
                        relationship_codes['Brother']['father'] = (member)
                    elif relationship == 'Mother':
                        relationship_codes['Proband']['mother'] = (member)
                        relationship_codes['Sister']['mother'] = (member)
                        relationship_codes['Brother']['mother'] = (member)

                    family_id = family
                    member_id = member
                    mother = relationship_codes[relationship]['mother']
                    father = relationship_codes[relationship]['father']
                    sex = sex_codes[clarity.get_sample(member).udf.get('Sex')]
                    phenotype = 0

                    line = '\t'.join([family_id, member_id, father, mother, sex, phenotype])
                    openfile.write(line)

        return ped_file

    def run_peddy(self, ped_file):
        cmd = 'peddy --plot --prefix %s %s %s' % (self.dataset, self.gatk_outfile, ped_file)

    def _run(self):
        self.run_gatk(self.gvcf_files)
        self.tabix_index()
        return self.run_gatk(self.gvcf_files) + self.tabix_index() + self.run_peddy(self.ped_file)

