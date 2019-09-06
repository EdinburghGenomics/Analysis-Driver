from egcg_core import executor
from analysis_driver import quality_control as qc
from analysis_driver.segmentation import Parameter
from analysis_driver.pipelines import common
from analysis_driver.pipelines.qc_gatk4 import SelectSNPs, SelectIndels, GATK4Stage, FastqIndex, SplitBWA, \
    MergeBamAndDup, SNPsFiltration, IndelsFiltration, MergeVariants, QCGATK4
from analysis_driver.pipelines.variant_calling_gatk4 import SplitHaplotypeCallerVC, \
    ScatterBaseRecalibrator, GatherBQSRReport, ScatterApplyBQSR, GatherRecalBam, GatherGVCF, VariantAnnotation, \
    GatherVCFVC, SplitGenotypeGVCFs


class VQSRFiltrationSNPs(GATK4Stage):
    """Run VariantRecalibrator on SNPs to create the truth and error models."""

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
    """Run VariantRecalibrator on Indels to create the truth and error models."""

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


class ApplyVQSR(GATK4Stage):
    """Apply truth and error models to the SNPs and INDELs"""

    input_vcf = Parameter()

    def _run(self):

        # Recalibrate the indels first
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
        recal_status = executor.execute(
            cmd,
            job_name='apply_vqsr_indels',
            working_dir=self.exec_dir,
            cpus=1,
            mem=16
        ).join()
        if recal_status == 0:
            # Then recalibrate the SNPs
            cmd = self.gatk_cmd(
                'ApplyVQSR',
                self.vqsr_filtered_vcf,
                memory=16,
                ext='-V {input_vcf} -mode SNP --tranches-file {ouput_tranches} --truth-sensitivity-filter-level 99.7 '
                    '--recal-file {recal_file}'.format(
                        input_vcf=self.vqsr_filtered_indels_vcf,
                        ouput_tranches=self.vqsr_snps_tranches,
                        recal_file=self.vqsr_snps_output_recall
                    )
            )
            recal_status = executor.execute(
                cmd,
                job_name='apply_vqsr_snps',
                working_dir=self.exec_dir,
                cpus=1,
                mem=16
            ).join()

        return recal_status


class HumanVarCallingGATK4(QCGATK4):
    _name = 'human_variant_calling_gatk4'

    def build(self):
        """Build the variant calling pipeline (for human)."""

        # merge fastq to do contamination check
        merge_fastqs = self.stage(common.MergeFastqs)
        contam = self.stage(qc.FastqScreen, previous_stages=[merge_fastqs], fq_pattern=merge_fastqs.fq_pattern)
        blast = self.stage(qc.Blast, previous_stages=[merge_fastqs], fastq_file=merge_fastqs.fq_pattern.replace('?', '1'))
        geno_val = self.stage(qc.GenotypeValidation, fq_pattern=merge_fastqs.fq_pattern, previous_stages=[merge_fastqs])

        # create fastq index then align and recalibrate via scatter-gather strategy
        fastq_index = self.stage(FastqIndex)
        split_bwa = self.stage(SplitBWA, previous_stages=[fastq_index])
        merge_bam_dup = self.stage(MergeBamAndDup, previous_stages=[split_bwa])
        base_recal = self.stage(ScatterBaseRecalibrator, previous_stages=[merge_bam_dup])
        gather_bqsr = self.stage(GatherBQSRReport, previous_stages=[base_recal])
        apply_bqsr = self.stage(ScatterApplyBQSR, previous_stages=[gather_bqsr])
        merge_bam = self.stage(GatherRecalBam, previous_stages=[apply_bqsr])

        # bam file QC
        verify_bam_id = self.stage(qc.VerifyBamID, bam_file=merge_bam.recal_bam, previous_stages=[merge_bam])
        samtools_stat = self.stage(common.SamtoolsStats, bam_file=merge_bam.recal_bam, previous_stages=[merge_bam])
        samtools_depth = self.stage(qc.SamtoolsDepth, bam_file=merge_bam.recal_bam, previous_stages=[merge_bam])

        # variants call via scatter-gather strategy
        haplotype_caller = self.stage(SplitHaplotypeCallerVC, bam_file=merge_bam.recal_bam, previous_stages=[merge_bam])
        gather_gcvf = self.stage(GatherGVCF, previous_stages=[haplotype_caller])
        genotype_gcvf = self.stage(SplitGenotypeGVCFs, previous_stages=[haplotype_caller])
        gather_vcf = self.stage(GatherVCFVC, previous_stages=[genotype_gcvf])

        variant_file = gather_vcf.genotyped_vcf
        steps_required = [gather_vcf]
        if 'snpEff' in self.dataset.genome_dict:
            # variant annotation
            annotate_vcf = self.stage(VariantAnnotation, previous_stages=[gather_vcf])
            variant_file = annotate_vcf.snps_effects_output_vcf
            steps_required = [annotate_vcf]

        # VQSR is disabled for now until we can make it perform better than hard filtering.

        # variant filtering with VQSR
        # filter_snps = self.stage(VQSRFiltrationSNPs, input_vcf=variant_file, previous_stages=steps_required)
        # filter_indels = self.stage(VQSRFiltrationIndels, input_vcf=variant_file, previous_stages=steps_required)
        # apply_vqsr = self.stage(ApplyVQSR, input_vcf=variant_file, previous_stages=[filter_indels, filter_snps])

        # variant filtering with Hard Filters
        select_snps = self.stage(SelectSNPs, input_vcf=variant_file, previous_stages=steps_required)
        select_indels = self.stage(SelectIndels, input_vcf=variant_file, previous_stages=steps_required)
        hard_filter_snps = self.stage(SNPsFiltration, previous_stages=[select_snps])
        hard_filter_indels = self.stage(IndelsFiltration, previous_stages=[select_indels])

        # Hard Filter variant merge
        merge_hf = self.stage(MergeVariants, stage_name='merge_variants_hard_filter',
                              vcf_files=[
                                  hard_filter_snps.hard_filtered_snps_vcf,
                                  hard_filter_indels.hard_filtered_indels_vcf
                              ],
                              output_vcf_file=hard_filter_indels.hard_filtered_vcf,
                              previous_stages=[hard_filter_snps, hard_filter_indels])

        sex_val = self.stage(qc.SexValidation, vcf_file=merge_hf.hard_filtered_vcf, previous_stages=[merge_hf])

        vcfstats = self.stage(qc.VCFStats, vcf_file=merge_hf.hard_filtered_vcf, previous_stages=[merge_hf])

        final_stages = [contam, blast, geno_val, sex_val, vcfstats, verify_bam_id, samtools_depth, samtools_stat,
                        gather_gcvf]

        output = self.stage(common.SampleDataOutput, previous_stages=final_stages, output_fileset='gatk4_human_var_calling')
        review = self.stage(common.SampleReview, previous_stages=[output])
        cleanup = self.stage(common.Cleanup, previous_stages=[review])

        return cleanup
