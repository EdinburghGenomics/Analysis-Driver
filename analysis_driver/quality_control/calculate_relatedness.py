import os
from collections import Counter
import re
import csv
from egcg_core import executor, util, clarity
from luigi import Parameter, ListParameter
from analysis_driver.config import default as cfg
from analysis_driver import segmentation
from analysis_driver.util.bash_commands import java_command
from analysis_driver.exceptions import PipelineError


class RelatednessStage(segmentation.Stage):
    @property
    def gatk_outfile(self):
        return os.path.join(self.job_dir, self.dataset.name + '_genotype_gvcfs.vcf')

    def family_id(self, sample_id):
        family_id = clarity.get_sample(sample_id).udf.get('Family ID')
        if not family_id:
            return 'No_ID'
        return family_id

    def relationship(self, member):
        relationship = clarity.get_sample(member).udf.get('Relationship')
        if not relationship:
            return 'Other'
        return relationship

    def sex(self, member):
        sex = clarity.get_sample(member).udf.get('Sex')
        if not sex:
            return 'No_Sex'
        return sex


class ParseRelatedness(RelatednessStage):
    parse_method = Parameter()
    ids = ListParameter()

    @property
    def user_sample_ids(self):
        user_to_internal_ids = {}
        for sample_id in self.ids:
            user_id = clarity.get_user_sample_name(sample_id)
            if user_id in user_to_internal_ids.keys():
                raise PipelineError('User ID %s appears more than once in sample list' % (user_id))
            user_to_internal_ids[user_id] = sample_id
        return user_to_internal_ids

    @property
    def reformat_relatedness_file(self):
        relatedness_reformat = (os.path.join(self.job_dir, self.dataset.name + '.relatedness2_reformat'))
        with open(os.path.join(self.job_dir, self.dataset.name + '.relatedness2')) as openfile:
            readfile = openfile.read()
            readfile_csv = re.sub('\t', ',', readfile)
            with open(relatedness_reformat, 'w') as outfile:
                outfile.write(readfile_csv)
        return relatedness_reformat

    @property
    def peddy_file(self):
        return os.path.join(self.job_dir, self.dataset.name + '.ped_check.csv')

    def write_all_fields(self, values):
        user_ids = self.user_sample_ids
        with open(os.path.join(self.job_dir, self.dataset.name + '.relatedness_output'), 'w') as outfile:
            for r in values:
                sample1 = r[0]
                sample2 = r[1]
                relatedness_values = r[2:]
                line = [sample1,
                        self.family_id(user_ids[sample1]),
                        self.relationship(user_ids[sample1]),
                        sample2,
                        self.family_id(user_ids[sample2]),
                        self.relationship(user_ids[sample1])]
                line.extend(relatedness_values)
                outfile.write('\t'.join(line) + '\n')

    def get_columns(self, filename, column_headers):
        columns_to_return = []
        with open(filename) as openfile:
            csvfile = csv.DictReader(openfile)
            for line in csvfile:
                columns_to_return.append([line[i] for i in column_headers])
        return columns_to_return

    def parse_all(self):
        exit_status = 0
        try:
            peddy_relatedness = self.get_columns(self.peddy_file, ['sample_a', 'sample_b', 'rel'])
            peddy_vcftools_relatedness = []
            with open(self.reformat_relatedness_file) as openfile:
                csvfile = csv.DictReader(openfile)
                for line in csvfile:
                    relatedness_samples = (Counter([line['INDV1'], line['INDV2']]))
                    [peddy_vcftools_relatedness.append([i[0],
                                                        i[1],
                                                        i[2],
                                                        line['RELATEDNESS_PHI']]) for i in peddy_relatedness if Counter([i[0], i[1]]) == relatedness_samples]
            self.write_all_fields(peddy_vcftools_relatedness)
        except (OSError, FileNotFoundError, NotADirectoryError) as e:
            self.error(str(e))
            exit_status += 1
        return exit_status

    def parse_single(self, sample1, sample2, relatedness, relatedness_file):
        exit_status = 0
        try:
            columns = self.get_columns(relatedness_file, [sample1, sample2, relatedness])
            self.write_all_fields(columns)
        except (OSError, FileNotFoundError, NotADirectoryError) as e:
            self.error(str(e))
            exit_status += 1
        return exit_status

    def _run(self):
        exit_status = 0
        if self.parse_method == 'parse_both':
            exit_status+=self.parse_all()
        elif self.parse_method == 'parse_relatedness':
            sample1 = 'INDV1'
            sample2 = 'INDV2'
            relatedness = 'RELATEDNESS_PHI'
            exit_status+=self.parse_single(sample1, sample2, relatedness, self.reformat_relatedness_file)
        elif self.parse_method == 'parse_peddy':
            sample1 = 'sample_a'
            sample2 = 'sample_b'
            relatedness = 'rel'
            exit_status+=self.parse_single(sample1, sample2, relatedness, self.peddy_file)
        return exit_status


class Genotype_gVCFs(RelatednessStage):
    gVCFs = ListParameter()
    reference = Parameter()

    def gatk_genotype_gvcfs_cmd(self):
        gvcf_variants = ' '. join(['--variant ' + util.find_file(i) for i in self.gVCFs])
        number_threads = 12
        return java_command(memory=50, tmp_dir=self.job_dir, jar=cfg['tools']['gatk']) + \
               ' -T GenotypeGVCFs -nt %s -R %s %s -o %s' % (
            number_threads, self.reference, gvcf_variants, self.gatk_outfile
        )

    def run_gatk(self):
        return executor.execute(
            self.gatk_genotype_gvcfs_cmd(),
            job_name='gatk_genotype_gvcfs',
            working_dir=self.job_dir,
            cpus=12,
            mem=50
        ).join()

    def _run(self):
        return self.run_gatk()


class Relatedness(RelatednessStage):

    def vcftools_relatedness_cmd(self):
        return '%s --relatedness2 --vcf %s --out %s' % (
            cfg['tools']['vcftools'], self.gatk_outfile, self.dataset.name
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


class Peddy(RelatednessStage):
    ids = ListParameter()

    @property
    def tabix_command(self):
         return "%s -f -p vcf %s" % (cfg['tools']['tabix'], self.gatk_outfile + '.gz')

    def bgzip(self):
       return executor.execute(
            self.bgzip_command,
            job_name='bgzip',
            working_dir=self.job_dir,
            cpus=1,
            mem=10
        ).join()

    @property
    def bgzip_command(self):
        return "%s %s" % (cfg['tools']['bgzip'], self.gatk_outfile)

    def tabix_index(self):
        return executor.execute(
            self.tabix_command,
            job_name='tabix',
            working_dir=self.job_dir,
            cpus=1,
            mem=10
        ).join()

    def write_ped_file(self):
        ped_file = os.path.join(self.job_dir, 'ped.fam')
        with open(ped_file, 'w') as openfile:
            for line in self.ped_file_content:
                openfile.write('\t'.join(line) + '\n')
        return ped_file

    def relationships(self, family):
        relationship_codes = {'Proband': {'Mother':'0', 'Father':'0'},
                              'Mother':{'Mother':'0', 'Father':'0'},
                              'Father':{'Mother':'0', 'Father':'0'},
                              'Sister':{'Mother':'0', 'Father':'0'},
                              'Brother': {'Mother':'0', 'Father':'0'},
                              'Other':{'Mother':'0', 'Father':'0'}}

        for member in family:
            relationship = self.relationship(member)
            if relationship in ['Father', 'Mother']:
                    for i in ['Proband', 'Sister', 'Brother']:
                        relationship_codes[i][relationship] = member
        return relationship_codes

    def get_member_details(self, family, all_families):
        family_lines = []
        family_info = self.relationships(all_families[family])
        for member in all_families[family]:
            sex_codes = {'Male': '1', 'Female': '2', 'No_Sex': '0'}
            relationship = self.relationship(member)
            member_id = clarity.get_user_sample_name(member)
            if not family_info[relationship]['Father'] == '0':
                father_id = clarity.get_user_sample_name(family_info[relationship]['Father'])
            else:
                father_id = family_info[relationship]['Father']
            if not family_info[relationship]['Mother'] == '0':
                mother_id = clarity.get_user_sample_name(family_info[relationship]['Mother'])
            else:
                mother_id = family_info[relationship]['Mother']

            line = [family,
                    member_id,
                   father_id,
                   mother_id,
                   sex_codes[self.sex(member)],
                   '0']
            family_lines.append(line)
        return family_lines

    def all_families(self):
        all_families = {}
        seen_user_samples = []
        for i in self.ids:
            if not clarity.get_user_sample_name(i) in seen_user_samples:
                seen_user_samples.append(clarity.get_user_sample_name(i))
                family_id = self.family_id(i)
                if not all_families.get(family_id):
                    all_families[family_id] = []
                all_families[family_id].append(i)
        return all_families

    @property
    def ped_file_content(self):
        ped_file_content = []
        all_families = self.all_families()
        for family in all_families:
            family_lines = self.get_member_details(family, all_families)
            for line in family_lines:
                ped_file_content.append(line)
        return ped_file_content

    @property
    def peddy_command(self):
        ped_file = self.write_ped_file()
        peddy_cmd = '%s --plot --prefix %s %s %s' % (cfg['tools']['peddy'], self.dataset.name, self.gatk_outfile + '.gz', ped_file)
        return peddy_cmd

    def run_peddy(self):
        return executor.execute(
            self.peddy_command,
            job_name='peddy',
            working_dir=self.job_dir,
            cpus=1,
            mem=10
        ).join()

    def _run(self):
        return self.bgzip() + self.tabix_index() + self.run_peddy()
