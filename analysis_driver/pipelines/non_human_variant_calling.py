import os
from analysis_driver import quality_control as qc, segmentation
from analysis_driver.util.bash_commands import java_command
from analysis_driver.tool_versioning import toolset
from analysis_driver.pipelines import common, variant_calling
from luigi import Parameter
from egcg_core import executor


class GATKStageNonHuman(segmentation.Stage):

    @property
    def gatk_run_dir(self):
        d = os.path.join(self.job_dir, 'gatk_var_calling')
        os.makedirs(d, exist_ok=True)
        return d

    def gatk_cmd(self, run_cls, output, input_bam=None, xmx=16, nct=16, ext=None):

        base_cmd = java_command(memory=xmx, tmp_dir=self.gatk_run_dir, jar=toolset['gatk']) + (
                    '-R {ref} -T {run_cls} --read_filter BadCigar --read_filter NotPrimaryAlignment '
                    '-o {output} -l INFO -U LENIENT_VCF_PROCESSING ')

        if input_bam:
            base_cmd += '-I {input_bam} '
        if ext:
            base_cmd += ext
        if nct > 1:
            base_cmd += ' -nct %s' % nct

        return base_cmd.format(
            ref=self.dataset.reference_genome,
            run_cls=run_cls,
            input_bam=input_bam,
            output=output
        )

    def select_snp_command(self):
        command = self.gatk_cmd('SelectVariants', self.raw_snp_vcf, nct=16)
        command += ' -V ' + self.raw_vcf
        command += ' -selectType SNP '
        return command

    def filter_snp_command(self):
        filter = [
            'QD < 2.0',
            'FS > 60.0',
            'MQ < 40.0',
            'MQRankSum < -12.5',
            'ReadPosRankSum < -8.0'
        ]
        filter = ' || '.join(filter)
        command = self.gatk_cmd('VariantFiltration', self.filter_snp_vcf, nct=16)
        command += ' -V ' + self.raw_snp_vcf
        command += ' --filterExpression ' + filter
        command += ' --filterName "SNP_FILTER"'
        return command

    def select_indel_command(self):
        command = self.gatk_cmd('SelectVariants', self.raw_indel_vcf, nct=16)
        command += ' -V ' + self.filter_snp_vcf
        command += ' -selectType INDEL '
        return command

    def filter_indel_vcf(self):
        filter = [
            'QD < 2.0',
            'FS > 200.0',
            'ReadPosRankSum < -20.0'
        ]
        filter = ' || '.join(filter)
        command = self.gatk_cmd('VariantFiltration', self.raw_vcf, nct=16)
        command += ' -V ' + self.raw_snp_vcf
        command += ' --filterExpression ' + filter
        command += ' --filterName "INDEL_FILTER"'
        return command

    @property
    def basename(self):
        return os.path.join(self.gatk_run_dir, self.dataset.user_sample_id)

    @property
    def sorted_bam(self):
        return os.path.join(self.job_dir, self.dataset.name + '.bam')

    @property
    def sample_gvcf(self):
        return self.basename + '.g.vcf'

    @property
    def raw_vcf(self):
        return self.basename + '.vcf'

    @property
    def raw_snp_vcf(self):
        return 'raw_snp_' + self.basename + '.vcf'

    @property
    def filter_snp_vcf(self):
        return 'filter_snp' + self.basename + '.vcf'

    @property
    def raw_indel_vcf(self):
        return 'raw_indel_' + self.basename + '.vcf'

    @property
    def filter_indel_vcf(self):
        return 'filter_indel_' + self.basename + '.vcf'

class HaplotypeCaller(GATKStageNonHuman):
    def _run(self):
        haplotype_cmd = self.gatk_cmd(
            'HaplotypeCaller',
            self.sample_gvcf,
            input_bam=self.sorted_bam,
            xmx=48,
            ext=('--pair_hmm_implementation VECTOR_LOGLESS_CACHING -ploidy 2 --emitRefConfidence GVCF '
                 '--variant_index_type LINEAR --variant_index_parameter 128000 --dbsnp ' + self.dbsnp)
        )
        for annot in ('BaseQualityRankSumTest', 'FisherStrand', 'GCContent', 'HaplotypeScore',
                      'HomopolymerRun', 'MappingQualityRankSumTest', 'MappingQualityZero', 'QualByDepth',
                      'ReadPosRankSumTest', 'RMSMappingQuality', 'DepthPerAlleleBySample', 'Coverage',
                      'ClippingRankSumTest', 'DepthPerSampleHC'):
            haplotype_cmd += ' --annotation ' + annot
        return executor.execute(
            haplotype_cmd,
            job_name='gatk_haplotype_call',
            working_dir=self.gatk_run_dir,
            cpus=16,
            mem=64
        ).join()


class GenotypeGVCFs(GATKStageNonHuman):
    def _run(self):
        genotype_gvcfs_cmd = self.gatk_cmd('GenotypeGVCFs', self.genotyped_vcf, nct=16, ext=' --variant ' + self.sample_gvcf)
        return executor.execute(
            genotype_gvcfs_cmd,
            job_name='gatk_genotype_gvcfs',
            working_dir=self.gatk_run_dir,
            mem=16
        ).join()


class SelectVariants(GATKStageNonHuman):
    command = Parameter()
    def _run(self):
        return executor.execute(
            self.command,
            job_name='var_filtration',
            working_dir=self.gatk_run_dir,
            mem=16
        ).join()


class VariantFiltration(GATKStageNonHuman):
    command = Parameter()
    def _run(self):
        return executor.execute(
            self.command,
            job_name='var_filtration',
            working_dir=self.gatk_run_dir,
            mem=16
        ).join()


def build_pipeline(dataset):
    def stage(cls, **params):
        return cls(dataset=dataset, **params)

    bam_file_production = common.build_bam_file_production(dataset)
    haplotype = stage(variant_calling.HaplotypeCaller, previous_stages=[bam_file_production])
    genotype = stage(GenotypeGVCFs, previous_stages=[haplotype])
    select_snp = stage(SelectVariants, command=GATKStageNonHuman.select_snp_command(), previous_stages=[genotype])
    filter_snp = stage(VariantFiltration, command=GATKStageNonHuman.filter_snp_command(), previous_stages=[select_snp])
    select_indel = stage(SelectVariants, command=GATKStageNonHuman.select_indel_command(), previous_stages=[filter_snp])
    filter_indel = stage(VariantFiltration, command=GATKStageNonHuman.filter_snp_command(), previous_stages=[select_indel])
    bgzip = stage(variant_calling.BGZip, previous_stages=[filter_indel])
    tabix = stage(variant_calling.Tabix, previous_stages=[bgzip])
    vcfstats = stage(qc.VCFStats, vcf_file=variant_calling.GATKStageHuman.genotyped_vcf, previous_stages=[tabix])
    output = stage(common.SampleDataOutput, previous_stages=[vcfstats], output_fileset='gatk_var_calling')
    _cleanup = stage(common.Cleanup, previous_stages=[output])
    review = stage(common.SampleReview, previous_stages=[_cleanup])
    return review

