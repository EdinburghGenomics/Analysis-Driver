import os
from analysis_driver import executor, util
from analysis_driver.external_data import clarity
from analysis_driver.exceptions import PipelineError
from analysis_driver.config import default as cfg
from analysis_driver.reader.mapping_stats_parsers import parse_genotype_concordance
from .quality_control_base import QualityControl


class GenotypeValidation(QualityControl):
    """
    This class will perform the Genotype validation steps. It subclasses Thread, allowing it to run in the
    background.
    """
    def __init__(self, dataset, working_dir, fastq_files, vcf_file=None, check_neighbour=False,
                 check_project=False, list_samples=None):
        """
        :param list[str] fastq_files: fastq files to run genotype validation on
        :param str vcf_file:
        :param check_neighbour:
        :param check_project:
        :param list_samples:
        """
        super().__init__(dataset, working_dir)
        self.sample_id = self.dataset.name
        self.fastq_files = fastq_files
        if vcf_file:
            self.seq_vcf_file = vcf_file
        else:
            self.seq_vcf_file = os.path.join(self.working_dir, self.sample_id + '_genotype_validation.vcf.gz')
        self.check_neighbour = check_neighbour
        self.check_project = check_project
        self.list_samples = list_samples
        self.sample2genotype_validation = {}

    @staticmethod
    def _bwa_aln(fastq_files, sample_name, expected_output_bam, reference):
        """
        Contruct a command that will perform the alignment and duplicate removal using bwa aln.
        :param list fastq_files: 1 or 2 fastq files
        :param str sample_name: the name of the sample that should be added in the read group
        :param str expected_output_bam: the name of the bam file that will be created
        :param str reference: the path to the reference file that will be used to align
        :rtype: str
        :return: A pipe-separated bash command
        """
        bwa_bin = cfg.query('tools', 'bwa', ret_default='bwa')
        if len(fastq_files) == 2:
            command_aln1 = '%s aln %s %s' % (bwa_bin, reference, fastq_files[0])
            command_aln2 = '%s aln %s %s' % (bwa_bin, reference, fastq_files[1])
            command_bwa = "%s sampe -r '@RG\\tID:1\\tSM:%s' %s <(%s) <(%s) %s %s" % (bwa_bin, sample_name,
                                                                                     reference, command_aln1,
                                                                                     command_aln2,
                                                                                     fastq_files[0],
                                                                                     fastq_files[1])
        elif len(fastq_files) == 1:
            command_aln1 = '%s aln %s %s' % (bwa_bin, reference, fastq_files[0])
            command_bwa = "%s samse -r '@RG\\tID:1\\tSM:%s' %s <(%s) %s" % (bwa_bin, sample_name, reference,
                                                                            command_aln1, fastq_files[0])

        else:
            raise PipelineError('Bad number of fastqs: ' + str(fastq_files))

        command_samblaster = '%s --removeDups' % (cfg.query('tools', 'samblaster', ret_default='samblaster'))
        command_samtools = '%s view -F 4 -Sb -' % (cfg.query('tools', 'samtools', ret_default='samtools'))
        command_sambamba = '%s sort -t 16 -o %s /dev/stdin' % (
            cfg.query('tools', 'sambamba', ret_default='sambamba'), expected_output_bam
        )

        return ' | '.join([command_bwa, command_samblaster, command_samtools, command_sambamba])

    def _bwa_alignment(self):
        """
        Run bwa alignment for all fastq files against a synthetic genome.
        :rtype: list
        :return list of bam file containing the reads aligned.
        """
        expected_bam = os.path.join(self.working_dir, self.sample_id + '_geno_val.bam')

        command = self._bwa_aln(
            self.fastq_files,
            self.sample_id,
            expected_bam,
            cfg.query('genotype-validation', 'reference')
        )
        
        self.dataset.start_stage('genotype_validation_bwa')
        bwa_executor = executor.execute(
            command,
            job_name='alignment_bwa',
            working_dir=self.working_dir,
            cpus=4,
            mem=8
        )
        exit_status = bwa_executor.join()
        self.dataset.end_stage('genotype_validation_bwa', exit_status)

        return expected_bam

    def _snp_calling(self, bam_file):
        """
        Call SNPs using GATK as defined in the config file.
        :param bam_file: The file containing all reads aligned to the synthetic genome.
        :rtype: str
        :return a vcf file that contains the variant for all samples.
        """
        gatk_command = [
            'java -Xmx4G -jar %s' % cfg.query('tools', 'gatk'),
            '-T UnifiedGenotyper',
            '-nt 4',
            '-R %s' % cfg.query('genotype-validation', 'reference'),
            ' --standard_min_confidence_threshold_for_calling 30.0',
            '--standard_min_confidence_threshold_for_emitting 0',
            '-out_mode EMIT_ALL_SITES',
            '-I %s' % bam_file,
            '-o %s' % self.seq_vcf_file
        ]

        self.dataset.start_stage('genotype_validation_gatk')
        gatk_executor = executor.execute(
            ' '.join(gatk_command),
            prelim_cmds=util.bash_commands.export_env_vars(),
            job_name='snpcall_gatk',
            working_dir=self.working_dir,
            cpus=4,
            mem=4
        )
        exit_status = gatk_executor.join()
        self.dataset.end_stage('genotype_validation_gatk', exit_status)
        return self.seq_vcf_file

    def _vcf_validation(self, sample2genotype):
        """
        Validate SNPs against genotype data found in the genotypes_repository
        :param sample2genotype:
        :rtype: list
        :return list of files containing the results of the validation.
        """
        list_commands = []
        sample2genotype_validation = {}
        # Make sure the file exists
        assert os.path.isfile(self.seq_vcf_file)

        if not os.path.isfile(self.seq_vcf_file+'.tbi'):
            self._index_vcf_gz(self.seq_vcf_file)
        for sample_name in sample2genotype:
            validation_results = os.path.join(self.working_dir, sample_name + '_genotype_validation.txt')
            sample2genotype_validation[sample_name] = validation_results
            gatk_command = ['java -Xmx4G -jar %s' % cfg.query('tools', 'gatk'),
                            '-T GenotypeConcordance',
                            '-eval:VCF %s ' % self.seq_vcf_file,
                            '-comp:VCF %s ' % sample2genotype.get(sample_name),
                            '-R %s' % cfg.query('genotype-validation', 'reference'),
                            ' > %s' % validation_results]
            list_commands.append(' '.join(gatk_command))

        self.dataset.start_stage('validation_genotype_concordance')
        genotype_concordance_executor = executor.execute(
            *list_commands,
            prelim_cmds=util.bash_commands.export_env_vars(),
            job_name='genotype_concordance',
            working_dir=self.working_dir,
            cpus=4,
            mem=8,
            log_commands=False
        )
        exit_status = genotype_concordance_executor.join()
        self.dataset.end_stage('validation_genotype_concordance', exit_status)
        return sample2genotype_validation

    def _rename_expected_genotype(self, sample2genotype):
        """This function assume only one sample in the header"""
        list_commands = []
        for sample_name in sample2genotype:
            genotype_vcf = sample2genotype.get(sample_name)
            cmd = "{bcftools} reheader -s <(echo {sn}) {genovcf} > {genovcf}.tmp; mv {genovcf}.tmp {genovcf}"
            list_commands.append(cmd.format(bcftools=cfg.query('tools', 'bcftools'), sn=self.sample_id, genovcf=genotype_vcf))

        exit_status = executor.execute(
            *list_commands,
            job_name='genotype_rename',
            working_dir=self.working_dir,
            cpus=4,
            mem=8,
            log_commands=False
        ).join()
        return exit_status

    def _index_vcf_gz(self, vcf_file):
        command = '{tabix} -p vcf {vcf}'.format(tabix=cfg.query('tools', 'tabix'), vcf=vcf_file)
        exit_status = executor.execute(
            command,
            job_name='index_vcf',
            working_dir=self.working_dir,
            cpus=1,
            mem=4
        ).join()
        return exit_status

    def _merge_validation_results(self, sample2genotype_validation):
        all_sample_lines = {}
        headers = None
        table_type = None
        for sample in sample2genotype_validation:
            if sample is self.sample_id:
                genotype_validation = sample2genotype_validation.get(sample)
                table_type, headers, lines = parse_genotype_concordance(genotype_validation)
                for line in lines:
                    sp_line = line.split()
                    all_sample_lines[sp_line[0]] = line
            else:
                genotype_validation = sample2genotype_validation.get(sample)
                table_type, headers, lines = parse_genotype_concordance(genotype_validation)
                sp_line = lines[1].split()
                sp_line[0] = sample
                all_sample_lines[sample] = '\t'.join(sp_line)
        with open(os.path.join(self.working_dir, self.sample_id + '_genotype_validation.txt'), 'w') as f:
            f.write(table_type + '\n')
            f.write('\t'.join(headers.split()) + '\n')
            for sample in sorted(all_sample_lines):
                f.write(all_sample_lines.get(sample) + '\n')

    def _genotype_validation(self):
        """
        Perform validation for each of the samples from a run
        :rtype: list
        :return list of file containing the results of the validation.
        """
        if not os.path.isfile(self.seq_vcf_file):
            # Generate self.seq_vcf_file
            bam_file = self._bwa_alignment()
            self._snp_calling(bam_file)

        # Compare against the sample or the plate
        samples_names = {self.sample_id}
        if self.check_neighbour:
            samples_names.update(clarity.get_samples_arrived_with(self.sample_id))
            samples_names.update(clarity.get_samples_genotyped_with(self.sample_id))
            samples_names.update(clarity.get_samples_sequenced_with(self.sample_id))
        if self.check_project:
            project_id = clarity.find_project_name_from_sample(self.sample_id)
            samples_names.update(clarity.get_sample_names_from_project(project_id))
        if self.list_samples:
            samples_names.update(self.list_samples)
        # Get the expected genotype from the LIMS
        sample2genotype = {}
        # TODO: Change to use batch calls which should be much faster
        for sample_name in samples_names:
            genotype_vcf = os.path.join(self.working_dir, sample_name + '_expected_genotype.vcf')
            genotype_vcf = clarity.get_sample_genotype(sample_name, genotype_vcf)
            if genotype_vcf:
                sample2genotype[sample_name] = genotype_vcf

        if sample2genotype:
            self._rename_expected_genotype(sample2genotype)
            self.sample2genotype_validation = self._vcf_validation(sample2genotype)
            if len(self.sample2genotype_validation) > 1:
                self._merge_validation_results(self.sample2genotype_validation)

    def run(self):
        try:
            self._genotype_validation()
        except Exception as e:
            self.exception = e
            self.exit_status = 8

    def join(self, timeout=None):
        super().join(timeout=timeout)
        if self.exception:
            raise self.exception
        return self.seq_vcf_file, self.sample2genotype_validation.get(self.sample_id)
