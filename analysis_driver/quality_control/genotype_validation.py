import os
from analysis_driver import writer

__author__ = 'tcezard'
from analysis_driver.config import default as cfg
validation_cfg = cfg.get('genotype-validation')

def bwa_alignment(fastq_files, sample_name, work_directory):
    """"""

    output_bam_file = os.path.join(work_directory, sample_name+'.bam')
    writer.commands.align_with_bwa_aln(fastq_files, sample_name, output_bam_file,
                                       validation_cfg.get('reference'))



def genotype_validation(fastq_files):
    bam_file = bwa_alignment(fastq_files, sample_name, work_directory)
    vcf_file = SNPs_calling(bam_file)
    validation_results = vcf_validation(vcf_file)