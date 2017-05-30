import os
from egcg_core import executor, util, clarity
from luigi import Parameter, ListParameter
from analysis_driver.config import default as cfg
from analysis_driver.segmentation import Stage
from analysis_driver.util.bash_commands import java_command

class Genotype_gVCFs(Stage):
    gVCFs = ListParameter()
    reference = Parameter()

    @property
    def gatk_outfile(self):
        return os.path.join(self.job_dir, self.dataset.name + '_genotype_gvcfs.vcf')

    def gatk_genotype_gvcfs_cmd(self):
        gvcf_variants = ' '. join(['--variant ' + util.find_file(i) for i in self.gvcf_files])
        number_threads = 12
        return java_command(memory=20, tmp_dir=self.job_dir, jar=cfg['tools']['gatk']) + \
               ' -T GenotypeGVCFs -nt %s -R %s %s -o %s' % (
            number_threads, self.reference, gvcf_variants, self.gatk_outfile
        )

    def run_gatk(self):
        return executor.execute(
            self.gatk_genotype_gvcfs_cmd(),
            job_name='gatk_genotype_gvcfs',
            working_dir=self.job_dir,
            cpus=30,
            mem=50
        ).join()

    def _run(self):
        return self.run_gatk()


class Relatedness(Stage):

    @property
    def gatk_outfile(self):
        return os.path.join(self.job_dir, self.dataset.name + '_genotype_gvcfs.vcf')

    project_id = Parameter()

    def vcftools_relatedness_cmd(self):
        return '%s --relatedness2 --vcf %s --out %s' % (
            cfg['tools']['vcftools'], self.gatk_outfile, self.project_id
        )

    def run_vcftools(self):
        return executor.execute(
            self.vcftools_relatedness_cmd(),
            job_name='vcftools_relatedness',
            working_dir=self.job_dir,
            cpus=1,
            mem=10
        ).join()

    def _run(self):
        return self.run_vcftools()


class Peddy(Stage):
    ids = Parameter()

    @property
    def gatk_outfile(self):
        return os.path.join(self.job_dir, self.dataset.name + '_genotype_gvcfs.vcf')

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

    def write_ped_file(self):
        ped_file = 'ped.fam'
        with open(ped_file) as openfile:
            openfile.write(self.ped_file_content)
        return ped_file

    def family_id(self, sample_id):
        family_id = clarity.get_sample(sample_id).udf.get('Family ID')
        if not family_id:
            return 'No_ID'
        return family_id

    def relationship(self, member):
        relationship = clarity.get_sample(member).udf.get('Relationship')
        if not relationship:
            return 'No_Relationship'
        return relationship

    def sex(self, member):
        sex = clarity.get_sample(member).udf.get('Sex')
        if not sex:
            return 'No_Sex'
        return sex

    @property
    def ped_file_content(self):
        sex_codes = {'Male': '1', 'Female': '2', 'No_Sex': '0'}
        all_families = {}
        for i in self.ids:
            family_id = self.family_id(i)
            if not all_families.get(family_id):
                all_families[family_id] = []
            all_families[family_id].append(i)

        ped_file_content = []
        for family in all_families:
            relationship_codes = {'Proband': {'mother':'0', 'father':'0'},
                                  'Mother':{'mother':'0', 'father':'0'},
                                  'Father':{'mother':'0', 'father':'0'},
                                  'Sister':{'mother':'0', 'father':'0'},
                                  'Brother': {'mother':'0', 'father':'0'},
                                  'Other':{'mother':'0', 'father':'0'}}

            for member in all_families[family]:
                relationship = self.relationship(member)
                if relationship == 'Father':
                    relationship_codes['Proband']['father'] = member
                    relationship_codes['Sister']['father'] = member
                    relationship_codes['Brother']['father'] = member
                elif relationship == 'Mother':
                    relationship_codes['Proband']['mother'] = member
                    relationship_codes['Sister']['mother'] = member
                    relationship_codes['Brother']['mother'] = member

                family_id = family
                member_id = member
                mother = relationship_codes[relationship]['mother']
                father = relationship_codes[relationship]['father']
                sex = sex_codes[self.sex(member)]
                phenotype = '0'
                line = [family_id, member_id, father, mother, sex, phenotype]
                ped_file_content.append(line)
        return ped_file_content

    @property
    def peddy_command(self):
        ped_file = self.write_ped_file()
        peddy_cmd = 'peddy --plot --prefix %s %s %s' % (self.dataset.name, self.genotyped_gvcfs, ped_file)
        return peddy_cmd

    def run_peddy(self):
        return executor.execute(
            self.peddy_command,
            job_name='peddy',
            working_dir=self.job_dir,
            cpus=10,
            mem=10
        ).join()

    def _run(self):
        return self.tabix_index() + self.run_peddy()
