import os
from egcg_core import executor
from analysis_driver import quality_control as qc
from analysis_driver.pipelines import Pipeline, common
from analysis_driver.pipelines.qc_gatk4 import SplitHaplotypeCaller, PostAlignmentScatter, GATK4Stage, \
    FastqIndex, SplitBWA, MergeBamAndDup, GatherVCF, MergeVariants, SelectSNPs, SelectIndels, SNPsFiltration, \
    IndelsFiltration
from analysis_driver.segmentation import Parameter
from analysis_driver.tool_versioning import toolset
from analysis_driver.util.bash_commands import picard_command, java_command


class PostAlignmentScatterVC(PostAlignmentScatter):
    """Generic class providing ability to split the genome in chromosomes."""

    chunk_handler = Parameter()

    def split_base_recal_grp(self, chunk):
        return os.path.join(self.split_file_dir, self.dataset.name + '_base_recal_grp_%s.grp' % chunk)

    def split_recal_bam(self, chunk):
        return os.path.join(self.split_file_dir, self.dataset.name + '_recal_%s.bam' % chunk)

    def gvcf_per_chunk(self, chunk):
        return os.path.join(self.split_file_dir, self.dataset.name + '_haplotype_caller_%s-%s-%s.g.vcf.gz' % chunk)

    def vcf_per_chunk(self, chunk):
        return os.path.join(self.split_file_dir, self.dataset.name + '_genotype_gvcf_%s-%s-%s.vcf.gz' % chunk)


class ScatterBaseRecalibrator(PostAlignmentScatterVC):
    """Run BaseRecalibrator on each Chromosome"""

    def base_recalibrator_cmd(self, chrom_names):
        return self.gatk_cmd(
            'BaseRecalibrator', self.split_base_recal_grp(chrom_names[0]), input=self.sorted_bam,
            memory=6, ext='--known-sites  ' + self.dbsnp + ' --intervals ' + ' --intervals '.join(chrom_names)
        )

    def _run(self):
        return executor.execute(
            *[self.base_recalibrator_cmd(chrom_names)
              for chrom_names in self.chunk_handler.split_genome_chromosomes()],
            job_name='gatk_base_recal',
            working_dir=self.exec_dir,
            cpus=1,
            mem=6
        ).join()


class GatherBQSRReport(PostAlignmentScatterVC):
    """Accumulate all reports created by BaseRecalibrator"""

    def _run(self):
        bqsr_reports_list = os.path.join(self.split_file_dir, self.dataset.name + '_bqsr_reports.list')
        with open(bqsr_reports_list, 'w') as open_file:
            for chrom_names in self.chunk_handler.split_genome_chromosomes():
                open_file.write(self.split_base_recal_grp(chrom_names[0]) + '\n')

        gather_bqsr_status = executor.execute(
            self.gatk_cmd('GatherBQSRReports', self.output_grp, input=bqsr_reports_list, memory=6, reference=None),
            job_name='gather_bqsr',
            working_dir=self.exec_dir,
            cpus=1,
            mem=6
        ).join()

        return gather_bqsr_status


class ScatterApplyBQSR(PostAlignmentScatterVC):
    """Run ApplyBQSR on each Chromosome"""

    def apply_bqsr_cmd(self, chrom_names):
        return self.gatk_cmd(
            'ApplyBQSR', self.split_recal_bam(chrom_names[0]), input=self.sorted_bam,
            memory=6, ext='--bqsr-recal-file ' + self.output_grp + ' --jdk-deflater --jdk-inflater '
                          '--static-quantized-quals 10 --static-quantized-quals 20 --static-quantized-quals 30 '
                          '--static-quantized-quals 40' + ' --intervals ' + ' --intervals '.join(chrom_names)
        )

    def _run(self):
        return executor.execute(
            *[self.apply_bqsr_cmd(chrom_names)
              for chrom_names in self.chunk_handler.split_genome_chromosomes(with_unmapped=True)],
            job_name='apply_bqsr',
            working_dir=self.exec_dir,
            cpus=1,
            mem=6
        ).join()


class GatherRecalBam(PostAlignmentScatterVC):
    """Merge all recalibrated bam file created by ApplyBQSR."""

    def _run(self):
        bam_file_list = os.path.join(self.split_file_dir, self.dataset.name + '_recal_bam.list')
        with open(bam_file_list, 'w') as open_file:
            for chrom_names in self.chunk_handler.split_genome_chromosomes(with_unmapped=True):
                open_file.write(self.split_recal_bam(chrom_names[0]) + '\n')

        gather_bam_status = executor.execute(
            picard_command('GatherBamFiles', input_file=bam_file_list, output_file=self.recal_bam,
                           tmp_dir=self.create_tmp_dir(), memory=6, assume_sorted=False,
                           picard_params={'CREATE_INDEX': 'true'}),
            job_name='gather_recal_bam',
            working_dir=self.exec_dir,
            cpus=1,
            mem=6
        ).join()

        return gather_bam_status


class SplitHaplotypeCallerVC(PostAlignmentScatterVC, SplitHaplotypeCaller):
    """
    Run HaplotypeCaller on each chunk of genomes to create a GVCF file.
    PostAlignmentScatterVC provides the file paths and SplitHaplotypeCaller the functions.
    """

    def haplotype_caller_cmd(self, chunk, region_file):
        haplotype_cmd = self.gatk_cmd(
            'HaplotypeCaller',
            self.gvcf_per_chunk(chunk),
            input=self.bam_file,
            memory=6,
            spark_core=1,
            ext=' --sample-ploidy 2 --emit-ref-confidence GVCF --intervals ' + region_file
        )
        if self.dbsnp:
            haplotype_cmd += ' --dbsnp ' + self.dbsnp
        if self.dataset.library_preparation == 'pcr-free':
            haplotype_cmd += '--pcr-indel-model NONE '

        return haplotype_cmd


class GatherGVCF(PostAlignmentScatterVC):
    """Collate all gvcf chunks into one."""

    def _run(self):
        gvcf_list = os.path.join(self.split_file_dir, self.dataset.name + '_g.vcf.list')
        with open(gvcf_list, 'w') as open_file:
            for chunks in self.chunk_handler.split_genome_in_chunks():
                open_file.write(self.gvcf_per_chunk(chunks[0]) + '\n')

        concat_vcf_status = executor.execute(
            self.gatk_picard_cmd('GatherVcfs', self.sample_gvcf, input=gvcf_list, memory=6, reference=None),
            job_name='gather_hap_call',
            working_dir=self.exec_dir,
            cpus=1,
            mem=6
        ).join()

        if concat_vcf_status == 0:
            concat_vcf_status = common.tabix_vcf(self.exec_dir, self.sample_gvcf)

        return concat_vcf_status


class SplitGenotypeGVCFs(PostAlignmentScatterVC):
    """Run GenotypeGVCFs on each chunk of the genome to create a VCF file."""

    chunk_handler = Parameter()

    def genotypegvcf_cmd(self, chunk, region_file):
        genotypegvcf_cmd = self.gatk_cmd(
            'GenotypeGVCFs',
            self.vcf_per_chunk(chunk),
            memory=6,
            spark_core=1,
            ext='--variant ' + self.gvcf_per_chunk(chunk) + ' --sample-ploidy 2 --intervals ' + region_file
        )
        if self.dbsnp:
            genotypegvcf_cmd += ' --dbsnp ' + self.dbsnp
        for annot in ('BaseQualityRankSumTest', 'ClippingRankSumTest', 'Coverage', 'DepthPerAlleleBySample',
                      'DepthPerSampleHC', 'FisherStrand', 'MappingQuality', 'MappingQualityRankSumTest',
                      'MappingQualityZero', 'QualByDepth', 'ReadPosRankSumTest', 'RMSMappingQuality'):
            genotypegvcf_cmd += ' --annotation ' + annot

        return genotypegvcf_cmd

    def _run(self):
        return executor.execute(
            *[self.genotypegvcf_cmd(chunk, region_file)
              for chunk, region_file in self.chunk_handler.split_genome_files(self.split_file_dir).items()],
            job_name='split_genotype_call',
            working_dir=self.exec_dir,
            cpus=1,
            mem=6
        ).join()


class GatherVCFVC(PostAlignmentScatterVC, GatherVCF):
    """
    Collate all vcf chunks into one.
    PostAlignmentScatterVC provide the file name and GatherVCF the functions.
    """
    pass


class VariantAnnotation(GATK4Stage):
    """Annotate a vcf file using snpEff."""

    def _run(self):
        cmd = 'set -o pipefail; ' + java_command(20, self.create_tmp_dir(), toolset['snpEff']) + \
              ('eff -hgvs -noLog -i vcf -o vcf '
               '-csvStats {effects_csv} -s {effects_html} {database_name} {input_vcf}'
               ' | {bgzip} --threads 16 -c > {output_vcf}').format(
                  effects_csv=self.snps_effects_csv,
                  effects_html=self.snps_effects_html,
                  database_name=self.dataset.genome_dict['snpEff'],  # to avoid getting None
                  input_vcf=self.genotyped_vcf,
                  bgzip=toolset['bgzip'],
                  output_vcf=self.snps_effects_output_vcf
              )
        snpeff_status = executor.execute(
            cmd,
            job_name='snpeff',
            working_dir=self.exec_dir,
            cpus=2,
            mem=20
        ).join()
        if snpeff_status == 0:
            snpeff_status = common.tabix_vcf(self.exec_dir, self.snps_effects_output_vcf)
        return snpeff_status


class VarCallingGATK4(Pipeline):
    toolset_type = 'gatk4_sample_processing'
    name = 'variant_calling_gatk4'

    def build(self):
        """Build the variant calling pipeline (for non human)."""
    
        # merge fastq to do contamination check
        merge_fastqs = self.stage(common.MergeFastqs)
        contam = self.stage(qc.FastqScreen, previous_stages=[merge_fastqs], fq_pattern=merge_fastqs.fq_pattern)
        blast = self.stage(qc.Blast, previous_stages=[merge_fastqs], fastq_file=merge_fastqs.fq_pattern.replace('?', '1'))
    
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
    
        # Variant stats
        vcfstats = self.stage(qc.VCFStats, vcf_file=merge_hf.hard_filtered_vcf, previous_stages=[merge_hf])
    
        final_stages = [contam, blast, vcfstats, verify_bam_id, samtools_depth, samtools_stat,
                        gather_gcvf, merge_hf]
    
        output = self.stage(common.SampleDataOutput, previous_stages=final_stages, output_fileset='gatk4_var_calling')
        review = self.stage(common.SampleReview, previous_stages=[output])
        cleanup = self.stage(common.Cleanup, previous_stages=[review])
    
        return cleanup
