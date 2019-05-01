from egcg_core import executor
from egcg_core.config import cfg

from analysis_driver import quality_control as qc
from analysis_driver.pipelines import common
from analysis_driver.pipelines.common import MergeFastqs, SamtoolsStats
from analysis_driver.pipelines.qc_gatk4 import SelectSNPs, SelectIndels, GATK4Stage, FastqIndex, SplitBWA, \
    MergeBamAndDup, SNPsFiltration, IndelsFiltration, MergeVariants
from analysis_driver.pipelines.variant_calling_gatk4 import SplitHaplotypeCallerVC, \
    ScatterBaseRecalibrator, GatherBQSRReport, ScatterApplyBQSR, GatherRecalBam, GatherGVCF, VariantAnnotation, \
    GatherVCFVC, SplitGenotypeGVCFs
from analysis_driver.segmentation import Parameter

toolset_type = 'gatk4_sample_processing'
name = 'human_variant_calling_gatk4'


class VQSRFiltrationSNPs(GATK4Stage):
    """Run VariantRecalibrator on SNPs to create the error model."""

    input_vcf = Parameter()

    def _run(self):
        vqsr_datasets = self.vqsr_datasets
        cmd = self.gatk_cmd(
            'VariantRecalibrator',
            self.vqsr_snps_output_recall,
            memory=16,
            ext='-V {input_vcf} '
                '--resource:hapmap,known=false,training=true,truth=true,prior=15.0 {hapmap} '
                '--resource:omni,known=false,training=true,truth=false,prior=12.0 {omni} '
                '--resource:1000G,known=false,training=true,truth=false,prior=10.0 {thousand_genomes} '
                '--resource:dbsnp,known=true,training=false,truth=false,prior=2.0 {dbsnp} '
                '--max-gaussians 4 -tranche 100.0 -tranche 99.99 -tranche 99.98 -tranche 99.97 -tranche 99.96 '
                '-tranche 99.95 -tranche 99.94 -tranche 99.93 -tranche 99.92 -tranche 99.91 -tranche 99.9 '
                '-tranche 99.8 -tranche 99.7 -tranche 99.6 -tranche 99.5 -tranche 99.0 -tranche 98.0 -tranche 90.0 '
                '-an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR -mode SNP '
                '--tranches-file {ouput_tranches} --rscript-file {output_R_script}'.format(
                    input_vcf=self.input_vcf,
                    hapmap=vqsr_datasets.get('hapmap'),
                    omni=vqsr_datasets.get('omni'),
                    thousand_genomes=vqsr_datasets.get('thousand_genomes'),
                    dbsnp=vqsr_datasets.get('dbsnp'),
                    ouput_tranches=self.vqsr_snps_tranches,
                    output_R_script=self.vqsr_snps_r_script
            )
        )
        return executor.execute(
            cmd,
            job_name='vqsr_snp',
            working_dir=self.exec_dir,
            cpus=1,
            mem=16
        ).join()


class VQSRFiltrationIndels(GATK4Stage):
    """Run VariantRecalibrator on Indels to create the error model."""

    input_vcf = Parameter()

    def _run(self):
        vqsr_datasets = self.vqsr_datasets
        cmd = self.gatk_cmd(
            'VariantRecalibrator',
            self.vqsr_indels_output_recall,
            memory=16,
            ext='-V {input_vcf} '
                '--resource:mills,known=false,training=true,truth=true,prior=12.0 {mills} '
                '--resource:dbsnp,known=true,training=false,truth=false,prior=2.0 {dbsnp} '
                '--max-gaussians 4 -tranche 100.0 -tranche 99.99 -tranche 99.98 -tranche 99.97 -tranche 99.96 '
                '-tranche 99.95 -tranche 99.94 -tranche 99.93 -tranche 99.92 -tranche 99.91 -tranche 99.9 '
                '-tranche 99.8 -tranche 99.7 -tranche 99.6 -tranche 99.5 -tranche 99.0 -tranche 98.0 -tranche 90.0 '
                '-an QD -an MQRankSum -an ReadPosRankSum -an FS -an SOR -an DP -mode INDEL '
                '--tranches-file {ouput_tranches} --rscript-file {output_R_script}'.format(
                    input_vcf=self.input_vcf,
                    mills=vqsr_datasets.get('mills'),
                    dbsnp=vqsr_datasets.get('dbsnp'),
                    ouput_tranches=self.vqsr_indels_tranches,
                    output_R_script=self.vqsr_indels_r_script,
                )
        )
        return executor.execute(
            cmd,
            job_name='vqsr_indel',
            working_dir=self.exec_dir,
            cpus=1,
            mem=16
        ).join()


class ApplyVQSRSNPs(GATK4Stage):
    """Apply error models to the SNPs."""

    input_vcf = Parameter()

    def _run(self):
        cmd = self.gatk_cmd(
            'ApplyVQSR',
            self.vqsr_filtered_snps_vcf,
            memory=16,
            ext='-V {input_vcf} -mode SNP --tranches-file {ouput_tranches} --truth-sensitivity-filter-level 99.0 '
                '--recal-file {recal_file}'.format(
                input_vcf=self.input_vcf,
                ouput_tranches=self.vqsr_snps_tranches,
                recal_file=self.vqsr_snps_output_recall
            )
        )
        return executor.execute(
            cmd,
            job_name='apply_vqsr_snps',
            working_dir=self.exec_dir,
            cpus=1,
            mem=16
        ).join()


class ApplyVQSRIndels(GATK4Stage):
    """Apply error models to the Indels."""

    input_vcf = Parameter()

    def _run(self):
        cmd = self.gatk_cmd(
            'ApplyVQSR',
            self.vqsr_filtered_indels_vcf,
            memory=8,
            ext='-V {input_vcf} -mode INDEL --tranches-file {ouput_tranches} --truth-sensitivity-filter-level 99.0 '
                '--recal-file {recal_file}'.format(
                    input_vcf=self.input_vcf,
                    ouput_tranches=self.vqsr_indels_tranches,
                    recal_file=self.vqsr_indels_output_recall
                )
        )
        return executor.execute(
            cmd,
            job_name='apply_vqsr_indels',
            working_dir=self.exec_dir,
            cpus=1,
            mem=16
        ).join()


def build_pipeline(dataset):
    """Build the variant calling pipeline (for human)."""

    def stage(cls, **params):
        return cls(dataset=dataset, **params)

    # merge fastq to do contamination check
    merge_fastqs = stage(MergeFastqs)
    contam = stage(qc.FastqScreen, previous_stages=[merge_fastqs], fq_pattern=merge_fastqs.fq_pattern)
    blast = stage(qc.Blast, previous_stages=[merge_fastqs], fastq_file=merge_fastqs.fq_pattern.replace('?', '1'))
    geno_val = stage(qc.GenotypeValidation, fq_pattern=merge_fastqs.fq_pattern, previous_stages=[merge_fastqs])

    # create fastq index then align and recalibrate via scatter-gather strategy
    fastq_index = stage(FastqIndex)
    split_bwa = stage(SplitBWA, previous_stages=[fastq_index])
    merge_bam_dup = stage(MergeBamAndDup, previous_stages=[split_bwa])
    base_recal = stage(ScatterBaseRecalibrator, previous_stages=[merge_bam_dup])
    gather_bqsr = stage(GatherBQSRReport, previous_stages=[base_recal])
    apply_bqsr = stage(ScatterApplyBQSR, previous_stages=[gather_bqsr])
    merge_bam = stage(GatherRecalBam, previous_stages=[apply_bqsr])

    # bam file QC
    verify_bam_id = stage(qc.VerifyBamID, bam_file=merge_bam.recal_bam, previous_stages=[merge_bam])
    samtools_stat = stage(SamtoolsStats, bam_file=merge_bam.recal_bam, previous_stages=[merge_bam])
    samtools_depth = stage(qc.SamtoolsDepth, bam_file=merge_bam.recal_bam, previous_stages=[merge_bam])

    # variants call via scatter-gather strategy
    haplotype_caller = stage(SplitHaplotypeCallerVC, bam_file=merge_bam.recal_bam, previous_stages=[merge_bam])
    gather_gcvf = stage(GatherGVCF, previous_stages=[haplotype_caller])
    genotype_gcvf = stage(SplitGenotypeGVCFs, previous_stages=[haplotype_caller])
    gather_vcf = stage(GatherVCFVC, previous_stages=[genotype_gcvf])

    variant_file = gather_vcf.genotyped_vcf
    steps_required = [gather_vcf]
    if 'snpEff' in cfg.query('genomes', dataset.genome_version):
        # variant annotation
        annotate_vcf = stage(VariantAnnotation, previous_stages=[gather_vcf])
        variant_file = annotate_vcf.snps_effects_output_vcf
        steps_required = [annotate_vcf]

    # variant filtering with VQSR
    filter_snps = stage(VQSRFiltrationSNPs, input_vcf=variant_file, previous_stages=steps_required)
    apply_vqsr_snps = stage(ApplyVQSRSNPs, input_vcf=variant_file, previous_stages=[filter_snps])
    filter_indels = stage(VQSRFiltrationIndels, input_vcf=variant_file, previous_stages=steps_required)
    apply_vqsr_indels = stage(ApplyVQSRIndels, input_vcf=variant_file, previous_stages=[filter_indels])

    # variant filtering with Hard Filters
    select_snps = stage(SelectSNPs, input_vcf=variant_file, previous_stages=steps_required)
    select_indels = stage(SelectIndels, input_vcf=variant_file, previous_stages=steps_required)
    hard_filter_snps = stage(SNPsFiltration, previous_stages=[select_snps])
    hard_filter_indels = stage(IndelsFiltration, previous_stages=[select_indels])

    # Hard Filter variant merge
    merge_hf = stage(MergeVariants, stage_name='merge_variants_hard_filter',
                     vcf_files=[hard_filter_snps.hard_filtered_snps_vcf, hard_filter_indels.hard_filtered_indels_vcf],
                     output_vcf_file=hard_filter_indels.hard_filtered_vcf,
                     previous_stages=[hard_filter_snps, hard_filter_indels])

    # Final variant merge
    merge_vqsr = stage(MergeVariants, stage_name='merge_variants_vqsr',
                       vcf_files=[filter_snps.vqsr_filtered_snps_vcf, filter_indels.vqsr_filtered_indels_vcf],
                       output_vcf_file=filter_indels.vqsr_filtered_vcf,
                       previous_stages=[apply_vqsr_snps, apply_vqsr_indels])

    # VQSR variant merge
    gender_val = stage(qc.GenderValidation, vcf_file=merge_vqsr.vqsr_filtered_vcf, previous_stages=[merge_vqsr])

    vcfstats = stage(qc.VCFStats, vcf_file=merge_vqsr.vqsr_filtered_vcf, previous_stages=[merge_vqsr])

    final_stages = [contam, blast, geno_val, gender_val, vcfstats, verify_bam_id, samtools_depth, samtools_stat,
                    gather_gcvf, merge_hf]

    output = stage(common.SampleDataOutput, previous_stages=final_stages, output_fileset='gatk4_human_var_calling')
    review = stage(common.SampleReview, previous_stages=[output])
    # cleanup = stage(common.Cleanup, previous_stages=[output])

    return [review]
