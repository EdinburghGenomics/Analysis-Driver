import os
from analysis_driver import quality_control as qc
from luigi import Parameter
from egcg_core import executor
from analysis_driver import segmentation
from analysis_driver.pipelines import common
from analysis_driver.config import default as cfg
from analysis_driver.pipelines.common import bgzip_and_tabix
from analysis_driver.util.bash_commands import java_command
from analysis_driver.tool_versioning import toolset
from analysis_driver.exceptions import AnalysisDriverError

toolset_type = 'non_human_sample_processing'
name = 'variant_calling'


class GATKStage(segmentation.Stage):
    @property
    def gatk_run_dir(self):
        d = os.path.join(self.job_dir, 'gatk_var_calling')
        os.makedirs(d, exist_ok=True)
        return d

    def gatk_cmd(self, run_cls, output, input_bam=None, xmx=16, nct=16, nt=16, ext=None):
        base_cmd = java_command(memory=xmx, tmp_dir=self.gatk_run_dir, jar=toolset['gatk']) + (
                    '-R {ref} -T {run_cls} --read_filter BadCigar --read_filter NotPrimaryAlignment '
                    '-o {output} -l INFO -U LENIENT_VCF_PROCESSING')

        if input_bam:
            base_cmd += ' -I {input_bam}'
        if ext:
            base_cmd += ext
        if nct > 1:
            base_cmd += ' -nct %s' % nct
        if nt > 1:
            base_cmd += ' -nt %s' % nt

        return base_cmd.format(
            ref=self.dataset.reference_genome,
            run_cls=run_cls,
            input_bam=input_bam,
            output=output
        )

    @property
    def basename(self):
        return os.path.join(self.gatk_run_dir, self.dataset.user_sample_id)

    @property
    def sorted_bam(self):
        return os.path.join(self.job_dir, self.dataset.name + '.bam')

    @property
    def recal_bam(self):
        return self.basename + '_recal.bam'

    @property
    def output_grp(self):
        return self.basename + '.grp'

    @property
    def output_intervals(self):
        return self.basename + '.intervals'

    @property
    def indel_realigned_bam(self):
        return self.basename + '_indel_realigned.bam'

    @property
    def sample_gvcf(self):
        return self.basename + '.g.vcf'

    @property
    def genotyped_vcf(self):
        return self.basename + '.vcf'

    @property
    def raw_snp_vcf(self):
        return self.basename + '_raw_snp.vcf'

    @property
    def filter_snp_vcf(self):
        return self.basename + '_filter_snp.vcf'

    @property
    def dbsnp(self):
        dbsnp = cfg.query('genomes', self.dataset.genome_version, 'dbsnp')
        if not dbsnp and self.dataset.pipeline.name == 'variant_calling':
            raise AnalysisDriverError('Could not find dbsnp file for %s' % self.dataset.name)
        return dbsnp

    @property
    def known_indels(self):
        return cfg.query('genomes', self.dataset.genome_version, 'known_indels')


class BaseRecal(GATKStage):
    def _run(self):
        return executor.execute(
            self.gatk_cmd(
                'BaseRecalibrator', self.output_grp, input_bam=self.sorted_bam,
                xmx=48, nt=1, ext=' --knownSites ' + self.dbsnp
            ),
            job_name='gatk_base_recal',
            working_dir=self.gatk_run_dir,
            cpus=16,
            mem=64
        ).join()


class PrintReads(GATKStage):
    def _run(self):
        return executor.execute(
            self.gatk_cmd(
                'PrintReads', self.recal_bam, input_bam=self.sorted_bam,
                xmx=48, nt=1, ext=' -BQSR ' + self.output_grp
            ),
            job_name='gatk_print_reads',
            working_dir=self.gatk_run_dir,
            cpus=16,
            mem=64
        ).join()


class RealignTarget(GATKStage):
    def _run(self):
        realign_target_cmd = self.gatk_cmd(
            'RealignerTargetCreator', self.output_intervals,
            input_bam=self.recal_bam, xmx=32, nct=1, nt=1
        )

        if self.known_indels:
            realign_target_cmd += ' --known ' + self.known_indels

        return executor.execute(
            realign_target_cmd,
            job_name='gatk_realign_target',
            working_dir=self.gatk_run_dir,
            cpus=1,
            mem=32
        ).join()


class Realign(GATKStage):
    def _run(self):
        realign_cmd = self.gatk_cmd(
            'IndelRealigner',
            self.indel_realigned_bam,
            input_bam=self.recal_bam,
            xmx=32,
            nct=1,
            nt=1,
            ext=' -targetIntervals ' + self.output_intervals
        )
        if self.known_indels:
            realign_cmd += ' --knownAlleles ' + self.known_indels
        return executor.execute(
            realign_cmd,
            job_name='gatk_indel_realign',
            working_dir=self.gatk_run_dir,
            cpus=1,
            mem=32
        ).join()


class HaplotypeCaller(GATKStage):
    input_bam = Parameter()

    def _run(self):
        haplotype_cmd = self.gatk_cmd(
            'HaplotypeCaller',
            self.sample_gvcf,
            input_bam=self.input_bam,
            xmx=48,
            nt=1,
            ext=(' --pair_hmm_implementation VECTOR_LOGLESS_CACHING -ploidy 2 --emitRefConfidence GVCF '
                 '--variant_index_type LINEAR --variant_index_parameter 128000 ')
        )
        for annot in ('BaseQualityRankSumTest', 'FisherStrand', 'GCContent', 'HaplotypeScore',
                      'HomopolymerRun', 'MappingQualityRankSumTest', 'MappingQualityZero', 'QualByDepth',
                      'ReadPosRankSumTest', 'RMSMappingQuality', 'DepthPerAlleleBySample', 'Coverage',
                      'ClippingRankSumTest', 'DepthPerSampleHC'):
            haplotype_cmd += ' --annotation ' + annot
        if self.dbsnp:
            haplotype_cmd += ' --dbsnp ' + self.dbsnp

        haplotype_status = executor.execute(
            haplotype_cmd,
            job_name='gatk_haplotype_call',
            working_dir=self.gatk_run_dir,
            cpus=16,
            mem=64
        ).join()

        return haplotype_status + bgzip_and_tabix(self.gatk_run_dir, self.sample_gvcf)


class GenotypeGVCFs(GATKStage):
    def _run(self):
        genotype_gvcfs_cmd = self.gatk_cmd('GenotypeGVCFs', self.genotyped_vcf, nct=1, ext=' --variant ' + self.sample_gvcf + '.gz')
        genotype_status = executor.execute(
            genotype_gvcfs_cmd,
            job_name='gatk_genotype_gvcfs',
            working_dir=self.gatk_run_dir,
            cpus=1,
            mem=16
        ).join()

        return genotype_status + bgzip_and_tabix(self.gatk_run_dir, self.genotyped_vcf)


class SelectVariants(GATKStage):
    def _run(self):
        select_var_command = self.gatk_cmd('SelectVariants', self.raw_snp_vcf, nct=1, nt=16)
        select_var_command += ' -V ' + self.genotyped_vcf + '.gz'
        select_var_command += ' -selectType SNP '
        select_variants_status = executor.execute(
            select_var_command,
            job_name='var_filtration',
            working_dir=self.gatk_run_dir,
            cpus=1,
            mem=16
        ).join()
        return select_variants_status + bgzip_and_tabix(self.gatk_run_dir, self.raw_snp_vcf)


class VariantFiltration(GATKStage):
    def _run(self):
        filter_array = [
            'QD < 2.0',
            'FS > 60.0',
            'MQ < 40.0',
            'MQRankSum < -12.5',
            'ReadPosRankSum < -8.0'
        ]
        filters = "'" + ' || '.join(filter_array) + "'"
        var_filter_command = self.gatk_cmd('VariantFiltration', self.filter_snp_vcf, nct=1, nt=1)
        var_filter_command += " -V " + self.raw_snp_vcf + '.gz'
        var_filter_command += " --filterExpression " + filters
        var_filter_command += " --filterName 'SNP_FILTER'"
        variant_filter_status = executor.execute(
            var_filter_command,
            job_name='var_filtration',
            working_dir=self.gatk_run_dir,
            cpus=1,
            mem=16
        ).join()
        return variant_filter_status + bgzip_and_tabix(self.gatk_run_dir, self.filter_snp_vcf)


def build_pipeline(dataset):

    def stage(cls, **params):
        return cls(dataset=dataset, **params)

    bam_file_production = common.build_bam_file_production(dataset)
    base_recal = stage(BaseRecal, previous_stages=bam_file_production)
    print_reads = stage(PrintReads, previous_stages=[base_recal])
    realign_target = stage(RealignTarget, previous_stages=[print_reads])
    realign = stage(Realign, previous_stages=[realign_target])
    haplotype = stage(HaplotypeCaller, input_bam=realign.indel_realigned_bam, previous_stages=[realign])
    genotype = stage(GenotypeGVCFs, previous_stages=[haplotype])
    select_snp = stage(SelectVariants, previous_stages=[genotype])
    filter_snp = stage(VariantFiltration, previous_stages=[select_snp])
    vcfstats = stage(qc.VCFStats, vcf_file=filter_snp.filter_snp_vcf + '.gz', previous_stages=[filter_snp])
    output = stage(common.SampleDataOutput, previous_stages=[vcfstats], output_fileset='gatk_var_calling')
    _cleanup = stage(common.Cleanup, previous_stages=[output])
    review = stage(common.SampleReview, previous_stages=[_cleanup])
    return review
