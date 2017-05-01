import os
from egcg_core import executor
from analysis_driver.pipelines.common import link_results_files, output_data, cleanup, get_dbsnp, get_known_indels
from analysis_driver.pipelines.common import build_bam_file_production
from analysis_driver.config import default as cfg
from analysis_driver.reader.version_reader import write_versions_to_yaml
from analysis_driver import segmentation


class GATKStage(segmentation.Stage):
    @property
    def gatk_run_dir(self):
        d = os.path.join(self.job_dir, 'gatk_var_calling')
        os.makedirs(d, exist_ok=True)
        return d

    def gatk_cmd(self, run_cls, input_bam, output, xmx=16, nct=16, ext=None):
        base_cmd = ('java -Xmx{xmx}m -XX:+UseSerialGC -Djava.io.tmpdir={tmpdir} -jar {gatk} -R {ref} '
                    '-I {input_bam} -T {run_cls} --read_filter BadCigar --read_filter NotPrimaryAlignment '
                    '-o {output} -l INFO -U LENIENT_VCF_PROCESSING ')

        if ext:
            base_cmd += ext
        if nct > 1:
            base_cmd += ' -nct %s' % nct

        return base_cmd.format(
            xmx=str(xmx * 1000),
            tmpdir=self.gatk_run_dir,
            gatk=cfg['tools']['gatk'],
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
    def dbsnp(self):
        return get_dbsnp(self.dataset.genome_version)

    @property
    def known_indels(self):
        return get_known_indels(self.dataset.genome_version)


class BaseRecal(GATKStage):
    def _run(self):
        return executor.execute(
            self.gatk_cmd('BaseRecalibrator', self.sorted_bam, self.output_grp, xmx=48, ext='--knownSites ' + self.dbsnp),
            job_name='gatk_base_recal',
            working_dir=self.gatk_run_dir,
            cpus=16,
            mem=64
        ).join()


class PrintReads(GATKStage):
    def _run(self):
        return executor.execute(
            self.gatk_cmd('PrintReads', self.sorted_bam, self.recal_bam, xmx=48, ext=' -BQSR ' + self.output_grp),
            job_name='gatk_print_reads',
            working_dir=self.gatk_run_dir,
            cpus=16,
            mem=64
        ).join()


class RealignTarget(GATKStage):
    def _run(self):
        realign_target_cmd = self.gatk_cmd('RealignerTargetCreator', self.recal_bam, self.output_intervals, nct=1)
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
            self.recal_bam,
            self.indel_realigned_bam,
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
            self.indel_realigned_bam,
            self.sample_gvcf,
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


class BGZip(GATKStage):
    def _run(self):
        return executor.execute(
            '%s %s' % (cfg['tools']['bgzip'], self.sample_gvcf),
            job_name='bgzip',
            working_dir=self.gatk_run_dir,
            cpus=1,
            mem=8
        ).join()


class Tabix(GATKStage):
    def _run(self):
        return executor.execute(
            '%s -p vcf %s' % (cfg['tools']['tabix'], self.sample_gvcf + '.gz'),
            job_name='tabix',
            working_dir=self.gatk_run_dir,
            cpus=1,
            mem=8
        ).join()


class Output(segmentation.Stage):
    def _run(self):
        dir_with_linked_files = link_results_files(self.dataset.name, self.job_dir, 'gatk_var_calling')
        write_versions_to_yaml(os.path.join(dir_with_linked_files, 'program_versions.yaml'))
        return output_data(self.dataset, self.job_dir, self.dataset.name, dir_with_linked_files)


class Cleanup(segmentation.Stage):
    def _run(self):
        return cleanup(self.dataset.name)


def build_pipeline(dataset):

    def stage(cls, **params):
        return cls(dataset=dataset, **params)

    bam_file_production = build_bam_file_production(dataset)

    base_recal = stage(BaseRecal, previous_stages=[bam_file_production])
    print_reads = stage(PrintReads, previous_stages=[base_recal])
    realign_target = stage(RealignTarget, previous_stages=[print_reads])
    realign = stage(Realign, previous_stages=[realign_target])
    haplotype = stage(HaplotypeCaller, previous_stages=[realign])
    bgzip = stage(BGZip, previous_stages=[haplotype])
    tabix = stage(Tabix, previous_stages=[bgzip])

    output = stage(Output, previous_stages=[tabix])
    _cleanup = stage(Cleanup, previous_stages=[output])
    return _cleanup
