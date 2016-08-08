from os import makedirs
from os.path import join
from egcg_core import executor
from analysis_driver.pipeline import Stage
from analysis_driver.config import default as cfg


class GATKStage(Stage):
    @property
    def gatk_run_dir(self):
        return join(cfg['jobs_dir'], self.dataset_name, 'gatk_var_calling')

    @property
    def sorted_bam(self):
        return join(self.job_dir, self.dataset_name + '.bam')

    @property
    def basename(self):
        return join(self.gatk_run_dir, self.dataset.user_sample_id)

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
        return cfg['references'][self.dataset.species]['dbsnp']

    @property
    def known_indels(self):
        return cfg['references'][self.dataset.species].get('known_indels')

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
            ref=cfg['references'][self.dataset.species]['fasta'],
            run_cls=run_cls,
            input_bam=input_bam,
            output=output
        )

    def _run(self):
        raise NotImplementedError


class Setup(GATKStage):
    def _run(self):
        makedirs(self.gatk_run_dir, exist_ok=True)
        return 0


class BaseRecal(GATKStage):
    previous_stages = Setup

    def _run(self):
        return executor.execute(
            self.gatk_cmd(
                'BaseRecalibrator',
                self.sorted_bam, self.output_grp, xmx=48, ext='--knownSites ' + self.dbsnp),
            job_name='gatk_base_recal',
            working_dir=self.gatk_run_dir,
            cpus=16,
            mem=64
        ).join()


class PrintReads(GATKStage):
    previous_stages = BaseRecal

    def _run(self):
        return executor.execute(
            self.gatk_cmd('PrintReads', self.sorted_bam, self.recal_bam, xmx=48, ext=' -BQSR ' + self.output_grp),
            job_name='gatk_print_reads',
            working_dir=self.gatk_run_dir,
            cpus=16,
            mem=64
        ).join()


class RealignTarget(GATKStage):
    previous_stages = PrintReads

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
    previous_stages = RealignTarget

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


class Haplotype(GATKStage):
    previous_stages = Realign

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
    previous_stages = Haplotype

    def _run(self):
        return executor.execute(
            '%s %s' % (cfg['tools']['bgzip'], self.sample_gvcf),
            job_name='bgzip',
            working_dir=self.gatk_run_dir,
            cpus=1,
            mem=8
        ).join()


class Tabix(GATKStage):
    previous_stages = BGZip

    def _run(self):
        return executor.execute(
            '%s -p vcf %s' % (cfg['tools']['tabix'], self.sample_gvcf + '.gz'),
            job_name='tabix',
            working_dir=self.gatk_run_dir,
            cpus=1,
            mem=8
        ).join()
