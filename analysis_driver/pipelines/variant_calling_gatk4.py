import os

from egcg_core import executor
from egcg_core.config import cfg

from analysis_driver import quality_control as qc
from analysis_driver.pipelines import common
from analysis_driver.pipelines.common import tabix_vcf, MergeFastqs, SamtoolsStats
from analysis_driver.pipelines.qc_gatk4 import SplitHaplotypeCaller, PostAlignmentScatter, GATK4FilePath, GATK4Stage, \
    FastqIndex, SplitBWA, MergeBamAndDup, GatherVCF, MergeVariants, SelectSNPs, SelectIndels, SNPsFiltration, \
    IndelsFiltration
from analysis_driver.tool_versioning import toolset
from analysis_driver.util.bash_commands import picard_command, java_command

toolset_type = 'gatk4_sample_processing'
name = 'variant_calling_gatk4'


class PostAlignmentScatterVC(PostAlignmentScatter):

    def split_genome_chromosomes(self, with_unmapped=False):
        """
        Split the genome per chromosomes aggregating smaller chromosomes to similar length as the longuest chromsome
        Code inspired from GATK best practice workflow:
        https://github.com/gatk-workflows/broad-prod-wgs-germline-snps-indels/blob/190945e358a6ee7a8c65bacd7b28c66527383376/PairedEndSingleSampleWf.wdl#L969

        :return: list of list of chromosome names
        """
        os.makedirs(self.split_file_dir, exist_ok=True)
        fai_file = self.dataset.reference_genome + '.fai'
        with open(fai_file) as open_file:
            sequence_tuple_list = []
            for line in open_file:
                sp_line = line.strip().split()
                sequence_tuple_list.append((sp_line[0], int(sp_line[1])))
            longest_sequence = sorted(sequence_tuple_list, key=lambda x: x[1], reverse=True)[0][1]
        chunks = []
        current_chunks = []
        chunks.append(current_chunks)
        temp_size = 0
        for sequence_tuple in sequence_tuple_list:
            if temp_size + sequence_tuple[1] <= longest_sequence:
                temp_size += sequence_tuple[1]
                current_chunks.append(sequence_tuple[0])
            else:
                current_chunks = []
                chunks.append(current_chunks)
                current_chunks.append(sequence_tuple[0])
                temp_size = sequence_tuple[1]
        # add the unmapped sequences as a separate line to ensure that they are recalibrated as well
        if with_unmapped:
            chunks.append(['unmapped'])
        return chunks

    def split_base_recal_grp(self, chunk):
        return os.path.join(self.split_file_dir, self.dataset.name + '_base_recal_grp_%s.grp' % chunk)

    def split_recal_bam(self, chunk):
        return os.path.join(self.split_file_dir, self.dataset.name + '_recal_%s.bam' % chunk)

    def gvcf_per_chunk(self, chunk):
        return os.path.join(self.split_file_dir, self.dataset.name + '_haplotype_caller_%s-%s-%s.g.vcf.gz' % chunk)

    def vcf_per_chunk(self, chunk):
        return os.path.join(self.split_file_dir, self.dataset.name + '_genotype_gvcf_%s-%s-%s.vcf.gz' % chunk)


class ScatterBaseRecalibrator(PostAlignmentScatterVC):

    def base_recalibrator_cmd(self, chrom_names):
        return self.gatk_cmd(
            'BaseRecalibrator', self.split_base_recal_grp(chrom_names[0]), input=self.sorted_bam,
            memory=6, ext='--known-sites  ' + self.dbsnp + ' --intervals ' + ' --intervals '.join(chrom_names)
        )

    def _run(self):
        return executor.execute(
            *[self.base_recalibrator_cmd(chrom_names)
              for chrom_names in self.split_genome_chromosomes()],
            job_name='gatk_base_recal',
            working_dir=self.exec_dir,
            cpus=1,
            mem=6
        ).join()


class GatherBQSRReport(PostAlignmentScatterVC):
    def _run(self):
        bqsr_reports_list = os.path.join(self.split_file_dir, self.dataset.name + '_bqsr_reports.list')
        with open(bqsr_reports_list, 'w') as open_file:
            for chrom_names in self.split_genome_chromosomes():
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
              for chrom_names in self.split_genome_chromosomes(with_unmapped=True)],
            job_name='apply_bqsr',
            working_dir=self.exec_dir,
            cpus=1,
            mem=6
        ).join()


class GatherRecalBam(PostAlignmentScatterVC):

    def _run(self):
        bam_file_list = os.path.join(self.split_file_dir, self.dataset.name + '_recal_bam.list')
        with open(bam_file_list, 'w') as open_file:
            for chrom_names in self.split_genome_chromosomes(with_unmapped=True):
                open_file.write(self.split_recal_bam(chrom_names[0]) + '\n')

        gather_bam_status = executor.execute(
            picard_command('GatherBamFiles', input_file=bam_file_list, output_file=self.recal_bam,
                           tmp_dir=self.gatk_run_dir, memory=6, assume_sorted=False,
                           picard_params={'CREATE_INDEX': 'true'}),
            job_name='gather_recal_bam',
            working_dir=self.exec_dir,
            cpus=1,
            mem=6
        ).join()

        return gather_bam_status


class SplitHaplotypeCallerVC(PostAlignmentScatterVC, SplitHaplotypeCaller):

    def haplotype_caller_cmd(self, chunk, region_file):
        haplotype_cmd = self.gatk_cmd(
            'HaplotypeCaller',
            self.gvcf_per_chunk(chunk),
            input=self.recal_bam,
            memory=6,
            spark_core=1,
            ext=' --sample-ploidy 2 --emit-ref-confidence GVCF --intervals ' + region_file
        )
        for annot in ('BaseQualityRankSumTest', 'ClippingRankSumTest', 'Coverage', 'DepthPerAlleleBySample',
                      'DepthPerSampleHC', 'FisherStrand', 'MappingQuality', 'MappingQualityRankSumTest',
                      'MappingQualityZero', 'QualByDepth', 'ReadPosRankSumTest', 'RMSMappingQuality'):
            haplotype_cmd += ' --annotation ' + annot
        if self.dbsnp:
            haplotype_cmd += ' --dbsnp ' + self.dbsnp
        if self.dataset.library_preparation == 'pcr-free':
            haplotype_cmd += '--pcr-indel-model NONE '

        return haplotype_cmd


class GatherGVCF(PostAlignmentScatterVC):
    def _run(self):
        gvcf_list = os.path.join(self.split_file_dir, self.dataset.name + '_g.vcf.list')
        with open(gvcf_list, 'w') as open_file:
            for chunks in self.split_genome_in_chunks():
                open_file.write(self.gvcf_per_chunk(chunks[0]) + '\n')

        concat_vcf_status = executor.execute(
            self.gatk_picard_cmd('GatherVcfs', self.sample_gvcf, input=gvcf_list, memory=6, reference=None),
            job_name='gather_hap_call',
            working_dir=self.exec_dir,
            cpus=1,
            mem=6
        ).join()

        if concat_vcf_status == 0:
            concat_vcf_status = tabix_vcf(self.exec_dir, self.sample_gvcf)

        return concat_vcf_status


class SplitGenotypeGVCFs(PostAlignmentScatterVC):

    def genotypegvcf_cmd(self, chunk, region_file):
        genotypegvcf_cmd = self.gatk_cmd(
            'GenotypeGVCFs',
            self.vcf_per_chunk(chunk),
            memory=6,
            spark_core=1,
            ext='--variant ' + self.gvcf_per_chunk(chunk) + ' --sample-ploidy 2 --intervals ' + region_file
        )
        for annot in ('BaseQualityRankSumTest', 'ClippingRankSumTest', 'Coverage', 'DepthPerAlleleBySample',
                      'DepthPerSampleHC', 'FisherStrand', 'MappingQuality', 'MappingQualityRankSumTest',
                      'MappingQualityZero', 'QualByDepth', 'ReadPosRankSumTest', 'RMSMappingQuality'):
            genotypegvcf_cmd += ' --annotation ' + annot

        return genotypegvcf_cmd

    def _run(self):
        return executor.execute(
            *[self.genotypegvcf_cmd(chunk, region_file)
              for chunk, region_file in self.split_genome_files().items()],
            job_name='split_genotype_call',
            working_dir=self.exec_dir,
            cpus=1,
            mem=6
        ).join()


class GatherVCFVC(PostAlignmentScatterVC, GatherVCF):
    pass


class VariantAnnotation(GATK4Stage):
    """Annotate a vcf file using snpEff."""

    def _run(self):
        cmd = 'set -o pipefail; ' + java_command(20, self.job_dir, toolset['snpEff']) + \
              (' eff -hgvs -noLog -i vcf -o vcf '
               '-csvStats {effects_csv} -s {effects_html} {database_name} {input_vcf}'
               ' | {bgzip} --threads 16 -c > {output_vcf}').format(
                  effects_csv=self.snps_effects_csv,
                  effects_html=self.snps_effects_html,
                  database_name=cfg.query('genomes', self.dataset.genome_version, 'snpEff'),
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
            snpeff_status = tabix_vcf(self.exec_dir, self.snps_effects_output_vcf)
        return snpeff_status


def build_pipeline(dataset):

    def stage(cls, **params):
        return cls(dataset=dataset, **params)

    # merge fastq to do contamination check
    merge_fastqs = stage(MergeFastqs)
    contam = stage(qc.FastqScreen, previous_stages=[merge_fastqs], fq_pattern=merge_fastqs.fq_pattern)
    blast = stage(qc.Blast, previous_stages=[merge_fastqs], fastq_file=merge_fastqs.fq_pattern.replace('?', '1'))

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

    # variant annotation
    annotate_vcf = stage(VariantAnnotation, previous_stages=[gather_vcf])

    # variant filtering with Hard Filters
    select_snps = stage(SelectSNPs, input_vcf=annotate_vcf.snps_effects_output_vcf, previous_stages=[annotate_vcf])
    select_indels = stage(SelectIndels, input_vcf=annotate_vcf.snps_effects_output_vcf, previous_stages=[annotate_vcf])
    hard_filter_snps = stage(SNPsFiltration, previous_stages=[select_snps])
    hard_filter_indels = stage(IndelsFiltration, previous_stages=[select_indels])

    # Hard Filter variant merge
    merge_hf = stage(MergeVariants, stage_name='merge_variants_hard_filter',
                     vcf_files=[hard_filter_snps.hard_filtered_snps_vcf, hard_filter_indels.hard_filtered_indels_vcf],
                     output_vcf_file=hard_filter_indels.hard_filtered_vcf,
                     previous_stages=[hard_filter_snps, hard_filter_indels])

    # Variant stats
    vcfstats = stage(qc.VCFStats, vcf_file=merge_hf.hard_filtered_vcf, previous_stages=[merge_hf])

    final_stages = [contam, blast, vcfstats, verify_bam_id, samtools_depth, samtools_stat,
                    gather_gcvf, merge_hf]

    output = stage(common.SampleDataOutput, previous_stages=final_stages, output_fileset='gatk4_var_calling')
    review = stage(common.SampleReview, previous_stages=[output])
    # cleanup = stage(common.Cleanup, previous_stages=[output])

    return [review]
