import os
from threading import Thread
from analysis_driver.app_logging import AppLogger
from analysis_driver import executor
from analysis_driver.clarity import get_genotype_information_from_lims, get_plate_id_and_well_from_lims, \
    get_sample_names_from_plate_from_lims
from analysis_driver.exceptions import AnalysisDriverError
from analysis_driver.notification import default as ntf
from analysis_driver.config import default as cfg


class GenotypeValidation(AppLogger, Thread):
    """
    This class will perform the Genotype validation steps. It subclasses Thread, allowing it to run in the
    background.
    """
    def __init__(self, fastqs_files, sample_id, vcf_file=None, check_plate=False):
        """
        :param dict[str, list[str]] sample_to_fastqs: a dict linking sample ids to their fastq files
        :param str run_id: the id of the run these sample were sequenced on.
        """
        self.fastqs_files = fastqs_files
        self.sample_id = sample_id
        self.work_directory = os.path.join(cfg['jobs_dir'], self.sample_id)
        self.validation_cfg = cfg.get('genotype-validation')
        if vcf_file:
            self.seq_vcf_file = vcf_file
        else:
            self.seq_vcf_file = os.path.join(self.work_directory, self.sample_id + '_genotype_validation.vcf.gz')
        self.check_plate = check_plate
        self.sample2genotype_validation
        self.exception = None
        Thread.__init__(self)

    def _bwa_aln(self, fastq_files, sample_name, expected_output_bam,  reference):
        """
        Contruct a command that will perform the alignment and duplicate removal using bwa aln.
        :param list fastq_files: 1 or 2 fastq files
        :param str sample_name: the name of the sample that should be added in the read group
        :param str expected_output_bam: the name of the bam file that will be created
        :param str reference: the path to the reference file that will be used to align
        :rtype: str
        :return: A pipe-separated bash command
        """
        bwa_bin = self.validation_cfg.get('bwa', 'bwa')
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
            raise AnalysisDriverError('Bad number of fastqs: ' + str(fastq_files))

        command_samblaster = '%s --removeDups' % (self.validation_cfg.get('samblaster', 'samblaster'))
        command_samtools = '%s view -F 4 -Sb -' % (self.validation_cfg.get('samtools', 'samtools'))
        command_sambamba = '%s sort -t 16 -o  %s /dev/stdin' % (
            self.validation_cfg.get('sambamba', 'sambamba'), expected_output_bam
        )

        return ' | '.join([command_bwa, command_samblaster, command_samtools, command_sambamba])

    def _bwa_alignment(self):
        """
        Run bwa alignment for all fastq files against a synthetic genome.
        :rtype: list
        :return list of bam file containing the reads aligned.
        """
        expected_bam = os.path.join(self.work_directory, self.sample_id + '_geno_val.bam')

        command = self._bwa_aln(
            self.fastqs_files,
            self.sample_id,
            expected_bam,
            self.validation_cfg.get('reference')
        )
        
        ntf.start_stage('genotype_validation_bwa')
        bwa_executor = executor.execute(
            [command],
            job_name='alignment_bwa',
            run_id=self.sample_id,
            cpus=4,
            mem=8
        )
        exit_status = bwa_executor.join()
        ntf.end_stage('genotype_validation_bwa', exit_status)

        return expected_bam

    def _snp_calling(self, bam_file):
        """
        Call SNPs using GATK as defined in the config file.
        :param bam_files: The file containing all the read aligned to the synthetic genome.
        :rtype: str
        :return a vcf file that contains the variant for all samples.
        """
        gatk_command = [
            'java -Xmx4G -jar %s' % self.validation_cfg.get('gatk'),
            '-T UnifiedGenotyper',
            '-nt 4',
            '-R %s' % self.validation_cfg.get('reference'),
            ' --standard_min_confidence_threshold_for_calling 30.0',
            '--standard_min_confidence_threshold_for_emitting 0',
            '-out_mode EMIT_ALL_SITES',
            '-I %s' % bam_file,
            '-o %s' % self.seq_vcf_file
        ]

        ntf.start_stage('genotype_validation_gatk')
        gatk_executor = executor.execute(
            [' '.join(gatk_command)],
            job_name='snpcall_gatk',
            run_id=self.sample_id,
            cpus=4,
            mem=4
        )
        exit_status = gatk_executor.join()
        ntf.end_stage('genotype_validation_gatk', exit_status)
        return self.seq_vcf_file

    def _vcf_validation(self, sample2genotype):
        """
        Validate SNPs against genotype data found in the genotypes_repository
        :param vcf_file: The vcf file containing the SNPs to validate.
        :rtype: list
        :return list of files containing the results of the validation.
        """
        list_commands = []
        sample2genotype_validation = {}
        for sample_name in sample2genotype:
            validation_results = os.path.join(self.work_directory, sample_name + '_genotype_validation.txt')
            sample2genotype_validation[sample_name] = validation_results
            gatk_command = ['java -Xmx4G -jar %s' % self.validation_cfg.get('gatk'),
                            '-T GenotypeConcordance',
                            '-eval:VCF %s ' % self.seq_vcf_file,
                            '-comp:VCF %s ' % sample2genotype.get(sample_name),
                            '-R %s' % self.validation_cfg.get('reference'),
                            ' > %s' % validation_results]
            list_commands.append(' '.join(gatk_command))


        ntf.start_stage('validation_genotype_concordance')
        genotype_concordance_executor = executor.execute(
            list_commands,
            job_name='genotype_concordance',
            run_id=self.sample_id,
            cpus=4,
            mem=8,
            log_command=False
        )
        exit_status = genotype_concordance_executor.join()
        ntf.end_stage('validation_genotype_concordance', exit_status)
        return sample2genotype_validation

    def _rename_expected_genotype(self, sample2genotype):
        """This function assume only one sample in the header"""
        list_commands = []
        for sample_name in sample2genotype:
            genotype_vcf = sample2genotype.get(sample_name)
            tmp_genotype = sample2genotype.get(sample_name) + '.tmp'
            cmd = "{bcftools} reheader -s <(echo {sn}) {genovcf} > {genovcf}.tmp; mv {genovcf}.tmp {genovcf}"
            list_commands.append(cmd.format(bcftools=self.validation_cfg.get('bcftools'), sn=sample_name, genovcf=genotype_vcf))

        exit_status = executor.execute(
            list_commands,
            job_name='genotype_rename',
            run_id=self.sample_id,
            cpus=4,
            mem=8,
            log_command=False
        ).join()
        return exit_status


    def _genotype_validation(self):
        """
        Perform validation for each of the samples from a run
        :rtype: list
        :return list of file containing the results of the validation.
        """
        if not os.path.isfile(self.seq_vcf_file):
            #Generate self.seq_vcf_file
            bam_file = self._bwa_alignment()
            self._snp_calling(bam_file)

        # Compare against the sample or the plate
        if self.check_plate:
            plate_id, well_id = get_plate_id_and_well_from_lims(self.sample_id)
            samples_names = get_sample_names_from_plate_from_lims(plate_id)
        else:
            samples_names = [self.sample_id]
        # Get the expected genotype from the LIMS
        sample2genotype = {}
        # TODO: Change to use batch calls which should be much faster
        for sample_name in samples_names:
            genotype_vcf = os.path.join(self.work_directory, sample_name + '_expected_genotype.vcf')
            genotype_vcf = get_genotype_information_from_lims(sample_name, genotype_vcf)
            if genotype_vcf:
                sample2genotype[sample_name] = genotype_vcf

        if sample2genotype:
            self._rename_expected_genotype(sample2genotype)
            self.sample2genotype_validation = self._vcf_validation(sample2genotype)

    def run(self):
        try:
            self._genotype_validation()
        except Exception as e:
            self.exception = e

    def join(self, timeout=None):
        super().join(timeout=timeout)
        if self.exception:
            raise self.exception
        return self.seq_vcf_file, self.sample2genotype_validation
