import os
from threading import Thread
from analysis_driver import writer
from analysis_driver.app_logging import AppLogger
from analysis_driver.executor import executor
from analysis_driver.notification import default as ntf
from analysis_driver.config import default as cfg


class GenotypeValidation(AppLogger, Thread):
    """ This class will perform the Genotype validation steps. It subclass Thread allowing it to be run in the background.
    """
    def __init__(self, sample_to_fastqs, run_id):
        """
        :param sample_to_fastqs: a dict linking the samples to their fastq file
        :param run_id the id of the run these sample were sequenced on.
        """
        self.sample_to_fastqs = sample_to_fastqs
        self.run_id = run_id
        self.validation_cfg = cfg.get('genotype-validation')
        Thread.__init__(self)

    def _align_with_bwa_aln(self, fastqs_files, sample_name, output_bam_file,  reference):
        """
        Contruct a command that will perform the alignment and duplicate removal using bwa aln.
        :param fastqs_files: array containing 1 or 2 fastq files
        :param sample_name: the name of the sample that should be added in the read group
        :param output_bam_file: the name of the bam file that will be created
        :param reference: the path to the reference file that will be used to align
        :rtype: str
        """
        bwa_bin = self.validation_cfg.get('bwa', 'bwa')
        if len(fastqs_files) == 2:
            command_aln1  = '%s aln %s %s' % (bwa_bin, reference, fastqs_files[0])
            command_aln2  = '%s aln %s %s' % (bwa_bin, reference, fastqs_files[1])
            command_bwa = "%s sampe -r '@RG\\tID:1\\tSM:%s' %s <(%s) <(%s) %s %s" % (bwa_bin, sample_name,
                                                                                     reference, command_aln1, command_aln2,
                                                                                     fastqs_files[0], fastqs_files[1])
        elif len(fastqs_files) == 1:
            command_aln1  = '%s aln %s %s' % (bwa_bin, reference, fastqs_files[0])
            command_bwa = "%s samse -r '@RG\\tID:1\\tSM:%s' %s <(%s) %s" % (bwa_bin, sample_name, reference,
                                                                            command_aln1, fastqs_files[0])

        command_samblaster = '%s --removeDups'%(self.validation_cfg.get('samblaster', 'samblaster'))
        command_samtools = '%s view -F 4 -Sb -'%(self.validation_cfg.get('samtools', 'samtools'))
        command_sambamba = '%s sort -t 16 -o  %s /dev/stdin'%(self.validation_cfg.get('sambamba', 'sambamba'), output_bam_file)

        return ' | '.join([command_bwa, command_samblaster, command_samtools, command_sambamba])


    def _bwa_alignment(self):
        """This function will run the alignment for all the fastq file against a syntetic genome.
        :rtype: list
        :return list of bam file containing the reads aligned."""
        list_commands = []
        list_output_bam = []
        work_dir = os.path.join(cfg['jobs_dir'], self.run_id)
        for sample_name in self.sample_to_fastqs:
            output_bam_file = os.path.join(work_dir, sample_name + '.bam')
            list_output_bam.append(output_bam_file)
            list_commands.append(self._align_with_bwa_aln(self.sample_to_fastqs.get(sample_name), sample_name,
                                                          output_bam_file, self.validation_cfg.get('reference')))
        bwa_writer = writer.get_script_writer('alignment_bwa', self.run_id, walltime=2, cpus=4, mem=8, jobs=len(list_commands))
        bwa_script = writer.write_jobs(
            bwa_writer,
            list_commands
        )

        ntf.start_stage('genotype_validation_bwa')
        bwa_executor = executor.ClusterExecutor(bwa_script, block=True)
        bwa_executor.start()
        exit_status = bwa_executor.join()
        ntf.end_stage('genotype_validation_bwa', exit_status)

        return list_output_bam

    def _SNPs_calling(self, bam_files):
        """
        This function will Call SNPs using GATK defined in the config file.
        :param bam_files: The file containing all the read aligned to the synthetic genome.
        :param run_id the id of the run these sample were sequenced on.
        :rtype: str
        :return a vcf file that contains the variant for all samples."""
        output_vcf_file = os.path.join(cfg['jobs_dir'], self.run_id, self.run_id + '_genotype_validation.vcf.gz')
        GATK_options = ['java -Xmx4G -jar %s' % self.validation_cfg.get('gatk'),
                        '-T UnifiedGenotyper',
                        '-t 4',
                        '-R %s' % self.validation_cfg.get('reference'),
                        ' --standard_min_confidence_threshold_for_calling 30.0',
                        '--standard_min_confidence_threshold_for_emitting 0',
                        '-out_mode EMIT_ALL_SITES']
        GATK_options.extend(['-I %s' % bam_file for bam_file in bam_files])
        GATK_options.append('-o %s' % output_vcf_file)

        gatk_writer = writer.get_script_writer('snpscall_gatk', self.run_id, walltime=2, cpus=4, mem=4, jobs=1)

        gatk_script = writer.write_jobs(
            gatk_writer,
            ' '.join(GATK_options)
        )
        ntf.start_stage('genotype_validation_gatk')
        gatk_executor = executor.ClusterExecutor(gatk_script, block=True)
        gatk_executor.start()
        exit_status = gatk_executor.join()
        ntf.end_stage('genotype_validation_gatk', exit_status)
        return output_vcf_file

    def _vcf_validation(self, vcf_file):
        """
        This will validate SNPs against genotype data found int the genotypes_repository
        :param vcf_file: The vcf file containing the SNPs to validate.
        :rtype: list
        :return list of file containing the results of the validation."""
        genotypes_dir = self.validation_cfg.get('genotypes_repository')
        work_dir = os.path.join(cfg['jobs_dir'], self.run_id)
        list_commands = []
        validation_results = []
        for sample_name in self.sample_to_fastqs:
            genotype_file = os.path.join(genotypes_dir, sample_name + '.vcf')
            validation_result = os.path.join(work_dir, sample_name + 'validation.txt')
            GATK_command = ['java -Xmx4G -jar %s' % self.validation_cfg.get('gatk'),
                            '-T GenotypeConcordance',
                            '-eval:VCF % ' % vcf_file,
                            '-comp:VCF % ' % genotype_file,
                            '-R %s' % self.validation_cfg.get('reference'),
                            ' > %s' % validation_result]
            list_commands.append(GATK_command)
            validation_results.append(validation_result)

        genotype_concordance_writer = writer.get_script_writer('validation_genotype_concordance', self.run_id,
                                                               walltime=2, cpus=4, mem=8, jobs=len(list_commands))
        genotype_concordance_writer = writer.write_jobs(genotype_concordance_writer,list_commands)

        ntf.start_stage('validation_genotype_concordance')
        bwa_executor = executor.ClusterExecutor(genotype_concordance_writer, block=True)
        bwa_executor.start()
        exit_status = bwa_executor.join()
        ntf.end_stage('validation_genotype_concordance', exit_status)
        return validation_results

    def _genotype_validation(self):
        """Perform the validation for each of the sample from a run
        :rtype: list
        :return list of file containing the results of the validation."""
        bam_files = self._bwa_alignment()
        vcf_file = self._SNPs_calling(bam_files)
        validation_results = self._vcf_validation(vcf_file)
        return validation_results

    def run(self):
        try :
            self.validation_results = self._genotype_validation()
        except Exception as e :
            self.exception = e

    def join(self, timeout=None):
        super().join(timeout=timeout)
        if self.exception:
            raise self.exception
        return self.validation_results
