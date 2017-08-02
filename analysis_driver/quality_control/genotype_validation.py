import os
from luigi import Parameter, BoolParameter, ListParameter
from egcg_core import executor, clarity, util
from analysis_driver.exceptions import PipelineError
from analysis_driver.util import bash_commands
from analysis_driver.config import default as cfg
from analysis_driver.reader.mapping_stats_parsers import parse_genotype_concordance
from analysis_driver.segmentation import Stage
from analysis_driver.tool_versioning import toolset


class GenotypeValidation(Stage):
    fq_pattern = Parameter()
    vcf_file = Parameter(default=None)
    check_neighbour = BoolParameter(default=False)
    check_project = BoolParameter(default=False)
    list_samples = ListParameter(default=[])
    validation_results = None

    @property
    def seq_vcf_file(self):
        return self.vcf_file or os.path.join(self.job_dir, self.dataset.name + '_genotype_validation.vcf.gz')

    @property
    def output_bam(self):
        """The Bam file that will be generated from BWA alignment"""
        return os.path.join(self.job_dir, self.dataset.name + '_geno_val.bam')

    def _gatk_command(self, memory):
        return bash_commands.java_command(memory=memory, tmp_dir=self.job_dir, jar=toolset['gatk'])

    def _bwa_aln(self, fastq_files):
        """
        Construct a command that will perform the alignment and duplicate removal using bwa aln.
        :param fastq_files: 1 or 2 fastqs
        :return: A pipe-separated bash command
        """
        header = '@RG\\tID:1\\tSM:' + self.dataset.name
        reference = cfg.query('genotype-validation', 'reference')

        if len(fastq_files) == 2:
            aln1 = '%s aln %s %s' % (toolset['bwa'], reference, fastq_files[0])
            aln2 = '%s aln %s %s' % (toolset['bwa'], reference, fastq_files[1])
            command_bwa = "%s sampe -r '%s' %s <(%s) <(%s) %s %s" % (
                toolset['bwa'], header, reference, aln1, aln2, fastq_files[0], fastq_files[1]
            )
        elif len(fastq_files) == 1:
            aln1 = '%s aln %s %s' % (toolset['bwa'], reference, fastq_files[0])
            command_bwa = "%s samse -r '%s' %s <(%s) %s" % (
                toolset['bwa'], header, reference, aln1, fastq_files[0]
            )

        else:
            raise PipelineError('Bad number of fastqs: ' + str(fastq_files))

        command_samblaster = '%s --removeDups' % toolset['samblaster']
        command_samtools = '%s view -F 4 -Sb -' % toolset['samtools']
        command_sambamba = '%s sort -t 16 -o %s /dev/stdin' % (
            toolset['sambamba'], self.output_bam
        )

        return ' | '.join([command_bwa, command_samblaster, command_samtools, command_sambamba])

    def _bwa_alignment(self):
        return executor.execute(
            self._bwa_aln(util.find_files(self.fq_pattern)),
            job_name='alignment_bwa',
            working_dir=self.job_dir,
            cpus=4,
            mem=8
        ).join()

    def _snp_calling(self, bam_file):
        """
        Call SNPs using GATK as defined in the config file.
        :param str bam_file: The file containing all reads aligned to the synthetic genome.
        """

        gatk_args = (
            '-T UnifiedGenotyper -nt 4 -R {ref} --standard_min_confidence_threshold_for_calling 30.0 '
            '--standard_min_confidence_threshold_for_emitting 0 -out_mode EMIT_ALL_SITES -I {bam_file} -o {vcf_file}'
        ).format(ref=cfg['genotype-validation']['reference'], bam_file=bam_file, vcf_file=self.seq_vcf_file)

        return executor.execute(
            self._gatk_command(memory=4) + gatk_args,
            prelim_cmds=bash_commands.export_env_vars(),
            job_name='snpcall_gatk',
            working_dir=self.job_dir,
            cpus=4,
            mem=4
        ).join()

    def _vcf_validation(self, sample2genotype):
        """
        Validate SNPs against genotype data found in the genotypes repository
        :param dict sample2genotype:
        :return files containing the results of the validation.
        """
        list_commands = []
        sample2genotype_validation = {}
        # Make sure the file exists and is indexed
        assert os.path.isfile(self.seq_vcf_file)
        if not os.path.isfile(self.seq_vcf_file + '.tbi'):
            self._index_vcf(self.seq_vcf_file)

        for sample_name in sample2genotype:
            validation_results = os.path.join(self.job_dir, sample_name + '_genotype_validation.txt')
            sample2genotype_validation[sample_name] = validation_results

            gatk_args = '-T GenotypeConcordance -eval:VCF {vcf} -comp:VCF {genotype} -R {ref} > {f}'.format(
                vcf=self.seq_vcf_file, genotype=sample2genotype.get(sample_name),
                ref=cfg['genotype-validation']['reference'], f=validation_results
            )
            list_commands.append(self._gatk_command(memory=4) + gatk_args)

        exit_status = executor.execute(
            *list_commands,
            prelim_cmds=bash_commands.export_env_vars(),
            job_name='genotype_concordance',
            working_dir=self.job_dir,
            cpus=4,
            mem=8,
            log_commands=False
        ).join()
        self.info('Checked vcfs for %s samples with exit status %s', len(sample2genotype), exit_status)
        return sample2genotype_validation

    def _rename_expected_genotype(self, sample2genotype):
        """This function assumes only one sample in the header"""
        base_cmd = '{bcftools} reheader -s <(echo {sn}) {genovcf} > {genovcf}.tmp; mv {genovcf}.tmp {genovcf}'
        commands = [
            base_cmd.format(bcftools=toolset['bcftools'], sn=self.dataset.name, genovcf=genovcf)
            for genovcf in sample2genotype.values()
        ]

        return executor.execute(
            *commands,
            job_name='genotype_rename',
            working_dir=self.job_dir,
            cpus=4,
            mem=8,
            log_commands=False
        ).join()

    def _index_vcf(self, vcf_file):
        return executor.execute(
            '{tabix} -p vcf {vcf}'.format(tabix=toolset['tabix'], vcf=vcf_file),
            job_name='index_vcf',
            working_dir=self.job_dir,
            cpus=1,
            mem=4
        ).join()

    def _merge_validation_results(self, sample2genotype_validation):
        all_sample_lines = {}
        headers = None
        table_type = None
        for sample, genotype_validation in sample2genotype_validation.items():
            table_type, headers, lines = parse_genotype_concordance(genotype_validation)
            if sample == self.dataset.name:
                for line in lines:
                    sp_line = line.split()
                    all_sample_lines[sp_line[0]] = line
            else:
                sp_line = lines[1].split()
                sp_line[0] = sample
                all_sample_lines[sample] = '\t'.join(sp_line)

        with open(os.path.join(self.job_dir, self.dataset.name + '_genotype_validation.txt'), 'w') as f:
            f.write(table_type + '\n')
            f.write('\t'.join(headers.split()) + '\n')
            for sample in sorted(all_sample_lines):
                f.write(all_sample_lines.get(sample) + '\n')

    def _run(self):
        if not os.path.isfile(self.seq_vcf_file):
            # Generate self.seq_vcf_file
            self._bwa_alignment()
            self._snp_calling(self.output_bam)

        # Compare against the sample or the plate
        sample_names = {self.dataset.name}
        if self.check_neighbour:
            sample_names.update(clarity.get_samples_arrived_with(self.dataset.name))
            sample_names.update(clarity.get_samples_genotyped_with(self.dataset.name))
            sample_names.update(clarity.get_samples_sequenced_with(self.dataset.name))
        if self.check_project:
            project_id = clarity.find_project_name_from_sample(self.dataset.name)
            sample_names.update(clarity.get_sample_names_from_project(project_id))
        if self.list_samples:
            sample_names.update(self.list_samples)

        # Get the expected genotype from the LIMS
        sample2genotype = {}
        # TODO: Change to use batch calls which should be much faster
        for s in sample_names:
            genotype_vcf = clarity.get_sample_genotype(
                s,
                os.path.join(self.job_dir, s + '_expected_genotype.vcf')
            )
            if genotype_vcf:
                sample2genotype[s] = genotype_vcf

        if sample2genotype:
            self._rename_expected_genotype(sample2genotype)
            sample2genotype_validation = self._vcf_validation(sample2genotype)
            if len(sample2genotype_validation) > 1:
                self._merge_validation_results(sample2genotype_validation)
            self.validation_results = sample2genotype_validation

        return 0
