import os
from analysis_driver import quality_control as qc
from egcg_core import executor
from analysis_driver import segmentation
from analysis_driver.pipelines import common
from analysis_driver.config import default as cfg
from analysis_driver.util.bash_commands import java_command
from analysis_driver.tool_versioning import toolset
from luigi import Parameter

toolset_type = 'non_human_sample_processing'
name = 'variant_calling'



class GATKStage(segmentation.Stage):
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

    @property
    def select_snp_command(self):
        command = self.gatk_cmd('SelectVariants', self.raw_snp_vcf, nct=16)
        command += ' -V ' + self.genotyped_vcf
        command += ' -selectType SNP '
        return command

    @property
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
    def raw_snp_vcf(self):
        return 'raw_snp_' + self.basename + '.vcf'

    @property
    def filter_snp_vcf(self):
        return 'filter_snp_' + self.basename + '.vcf'

    @property
    def dbsnp(self):
        return cfg.query('genomes', self.dataset.genome_version, 'dbsnp')

    @property
    def known_indels(self):
        return cfg.query('genomes', self.dataset.genome_version, 'known_indels')


class BaseRecal(GATKStage):
    def _run(self):
        return executor.execute(
            self.gatk_cmd('BaseRecalibrator', self.output_grp, input_bam=self.sorted_bam, xmx=48, ext='--knownSites ' + self.dbsnp),
            job_name='gatk_base_recal',
            working_dir=self.gatk_run_dir,
            cpus=16,
            mem=64
        ).join()


class PrintReads(GATKStage):
    def _run(self):
        return executor.execute(
            self.gatk_cmd('PrintReads', self.recal_bam, input_bam=self.sorted_bam, xmx=48, ext=' -BQSR ' + self.output_grp),
            job_name='gatk_print_reads',
            working_dir=self.gatk_run_dir,
            cpus=16,
            mem=64
        ).join()


class RealignTarget(GATKStage):
    def _run(self):
        realign_target_cmd = self.gatk_cmd('RealignerTargetCreator', self.output_intervals, input_bam=self.recal_bam, nct=1)
        if self.known_indels:
            realign_target_cmd += ' --known ' + self.known_indels

        return executor.execute(
            realign_target_cmd,
            job_name='gatk_realign_target',
            working_dir=self.gatk_run_dir,
            mem=16
        ).join()


class Realign(GATKStage):
    def _run(self):
        realign_cmd = self.gatk_cmd(
            'IndelRealigner',
            self.indel_realigned_bam,
            input_bam=self.recal_bam,
            nct=1,
            ext='-targetIntervals ' + self.output_intervals
        )
        if self.known_indels:
            realign_cmd += ' --knownAlleles ' + self.known_indels
        return executor.execute(
            realign_cmd,
            job_name='gatk_indel_realign',
            working_dir=self.gatk_run_dir,
            mem=16
        ).join()


class HaplotypeCaller(GATKStage):
    def _run(self):
        haplotype_cmd = self.gatk_cmd(
            'HaplotypeCaller',
            self.sample_gvcf,
            input_bam=self.indel_realigned_bam,
            xmx=48,
            ext=('--pair_hmm_implementation VECTOR_LOGLESS_CACHING -ploidy 2 --emitRefConfidence GVCF '
                 '--variant_index_type LINEAR --variant_index_parameter 128000 ')
        )
        for annot in ('BaseQualityRankSumTest', 'FisherStrand', 'GCContent', 'HaplotypeScore',
                      'HomopolymerRun', 'MappingQualityRankSumTest', 'MappingQualityZero', 'QualByDepth',
                      'ReadPosRankSumTest', 'RMSMappingQuality', 'DepthPerAlleleBySample', 'Coverage',
                      'ClippingRankSumTest', 'DepthPerSampleHC'):
            haplotype_cmd += ' --annotation ' + annot
            if self.dbsnp:
                haplotype_cmd += ' --dbsnp ' + self.dbsnp

        return executor.execute(
            haplotype_cmd,
            job_name='gatk_haplotype_call',
            working_dir=self.gatk_run_dir,
            cpus=16,
            mem=64
        ).join()


class GenotypeGVCFs(GATKStage):
    def _run(self):
        genotype_gvcfs_cmd = self.gatk_cmd('GenotypeGVCFs', self.genotyped_vcf, nct=16, ext=' --variant ' + self.sample_gvcf)
        return executor.execute(
            genotype_gvcfs_cmd,
            job_name='gatk_genotype_gvcfs',
            working_dir=self.gatk_run_dir,
            mem=16
        ).join()



class SelectVariants(GATKStage):
    command = Parameter()
    def _run(self):
        return executor.execute(
            self.command,
            job_name='var_filtration',
            working_dir=self.gatk_run_dir,
            mem=16
        ).join()


class VariantFiltration(GATKStage):
    command = Parameter()
    def _run(self):
        return executor.execute(
            self.command,
            job_name='var_filtration',
            working_dir=self.gatk_run_dir,
            mem=16
        ).join()


class BGZip(GATKStage):
    def _run(self):
        return executor.execute(
            '%s %s' % (toolset['bgzip'], self.filter_snp_vcf),
            job_name='bgzip',
            working_dir=self.gatk_run_dir,
            cpus=1,
            mem=8
        ).join()


class Tabix(GATKStage):
    def _run(self):
        return executor.execute(
            '%s -p vcf %s' % (toolset['tabix'], self.filter_snp_vcf + '.gz'),
            job_name='tabix',
            working_dir=self.gatk_run_dir,
            cpus=1,
            mem=8
        ).join()


def build_pipeline(dataset):

    def stage(cls, **params):
        return cls(dataset=dataset, **params)

    bam_file_production = common.build_bam_file_production(dataset)
    base_recal = stage(BaseRecal, previous_stages=bam_file_production)
    print_reads = stage(PrintReads, previous_stages=[base_recal])
    realign_target = stage(RealignTarget, previous_stages=[print_reads])
    realign = stage(Realign, previous_stages=[realign_target])
    haplotype = stage(HaplotypeCaller, previous_stages=[realign])
    genotype = stage(GenotypeGVCFs, previous_stages=[haplotype])
    select_snp = stage(SelectVariants, command=GATKStage.select_snp_command, previous_stages=[genotype])
    filter_snp = stage(VariantFiltration, command=GATKStage.filter_snp_command, previous_stages=[select_snp])
    vcfstats = stage(qc.VCFStats, vcf_file=GATKStage.filter_snp_vcf, previous_stages=[filter_snp])
    bgzip = stage(BGZip, previous_stages=[vcfstats])
    tabix = stage(Tabix, previous_stages=[bgzip])
    output = stage(common.SampleDataOutput, previous_stages=[tabix], output_fileset='gatk_var_calling')
    _cleanup = stage(common.Cleanup, previous_stages=[output])
    review = stage(common.SampleReview, previous_stages=[_cleanup])
    return review
