import os
from threading import Thread
from analysis_driver.app_logging import AppLogger
from analysis_driver import executor
from analysis_driver.exceptions import AnalysisDriverError
from analysis_driver.notification import default as ntf
from analysis_driver.config import default as cfg


class GenotypeValidation(AppLogger, Thread):
    """
    This class will perform the Genotype validation steps. It subclasses Thread, allowing it to run in the
    background.
    """
    def __init__(self, sample_to_fastqs, run_id):
        """
        :param dict[str, list[str]] sample_to_fastqs: a dict linking sample ids to their fastq files
        :param str run_id: the id of the run these sample were sequenced on.
        """
        self.sample_to_fastqs = sample_to_fastqs
        self.run_id = run_id
        self.validation_cfg = cfg.get('genotype-validation')
        self.validation_results = None
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
            raise AnalysisDriverError('Bad number of fastqs: ' + str(fastq_files))

        command_samblaster = '%s --removeDups' % (cfg.query('tools', 'samblaster', ret_default='samblaster'))
        command_samtools = '%s view -F 4 -Sb -' % (cfg.query('tools', 'samtools', ret_default='samtools'))
        command_sambamba = '%s sort -t 16 -o  %s /dev/stdin' % (
            cfg.query('tools', 'sambamba', ret_default='sambamba'), expected_output_bam
        )

        return ' | '.join([command_bwa, command_samblaster, command_samtools, command_sambamba])

    def _bwa_alignment(self):
        """
        Run bwa alignment for all fastq files against a synthetic genome.
        :rtype: list
        :return list of bam file containing the reads aligned.
        """
        commands = []
        output_bams = []
        work_dir = os.path.join(cfg['jobs_dir'], self.run_id)
        for sample_name in self.sample_to_fastqs:
            expected_bam = os.path.join(work_dir, sample_name + '.bam')
            output_bams.append(expected_bam)
            commands.append(
                self._bwa_aln(
                    self.sample_to_fastqs.get(sample_name),
                    sample_name,
                    expected_bam,
                    self.validation_cfg.get('reference')
                )
            )
        
        ntf.start_stage('genotype_validation_bwa')
        bwa_executor = executor.execute(
            commands,
            job_name='alignment_bwa',
            run_id=self.run_id,
            walltime=6,
            cpus=4,
            mem=8
        )
        exit_status = bwa_executor.join()
        ntf.end_stage('genotype_validation_bwa', exit_status)

        return output_bams

    def _snp_calling(self, bam_files):
        """
        Call SNPs using GATK as defined in the config file.
        :param bam_files: The file containing all the read aligned to the synthetic genome.
        :rtype: str
        :return a vcf file that contains the variant for all samples.
        """
        output_vcf = os.path.join(cfg['jobs_dir'], self.run_id, self.run_id + '_genotype_validation.vcf.gz')
        gatk_command = ['java -Xmx4G -jar %s' % cfg['tools']['gatk'],
                        '-T UnifiedGenotyper',
                        '-nt 4',
                        '-R %s' % self.validation_cfg.get('reference'),
                        ' --standard_min_confidence_threshold_for_calling 30.0',
                        '--standard_min_confidence_threshold_for_emitting 0',
                        '-out_mode EMIT_ALL_SITES']
        gatk_command.extend(['-I %s' % bam_file for bam_file in bam_files])
        gatk_command.append('-o %s' % output_vcf)

        ntf.start_stage('genotype_validation_gatk')
        gatk_executor = executor.execute(
            [' '.join(gatk_command)],
            job_name='snpcall_gatk',
            run_id=self.run_id,
            walltime=2,
            cpus=4,
            mem=4
        )
        exit_status = gatk_executor.join()
        ntf.end_stage('genotype_validation_gatk', exit_status)
        return output_vcf

    def _vcf_validation(self, vcf_file):
        """
        Validate SNPs against genotype data found in the genotypes_repository
        :param vcf_file: The vcf file containing the SNPs to validate.
        :rtype: list
        :return list of files containing the results of the validation.
        """
        genotypes_dir = self.validation_cfg.get('genotypes_repository')
        work_dir = os.path.join(cfg['jobs_dir'], self.run_id)
        list_commands = []
        validation_results = []
        for sample_name in self.sample_to_fastqs:
            genotype_file = os.path.join(genotypes_dir, sample_name + '.vcf')
            validation_result = os.path.join(work_dir, sample_name + 'validation.txt')
            gatk_command = ['java -Xmx4G -jar %s' % cfg['tools']['gatk'],
                            '-T GenotypeConcordance',
                            '-eval:VCF %s ' % vcf_file,
                            '-comp:VCF %s ' % genotype_file,
                            '-R %s' % self.validation_cfg.get('reference'),
                            ' > %s' % validation_result]
            list_commands.append(' '.join(gatk_command))
            validation_results.append(validation_result)

        ntf.start_stage('validation_genotype_concordance')
        genotype_concordance_executor = executor.execute(
            list_commands,
            job_name='genotype_concordance',
            run_id=self.run_id,
            walltime=2,
            cpus=4,
            mem=8
        )
        exit_status = genotype_concordance_executor.join()
        ntf.end_stage('validation_genotype_concordance', exit_status)
        return validation_results

    def _genotype_validation(self):
        """
        Perform validation for each of the samples from a run
        :rtype: list
        :return list of file containing the results of the validation.
        """
        bam_files = self._bwa_alignment()
        vcf_file = self._snp_calling(bam_files)
        validation_results = self._vcf_validation(vcf_file)
        return validation_results

    def run(self):
        try:
            self.validation_results = self._genotype_validation()
        except Exception as e:
            self.exception = e

    def join(self, timeout=None):
        super().join(timeout=timeout)
        if self.exception:
            raise self.exception
        return self.validation_results
