import os
from tests.test_analysisdriver import NamedMock, TestAnalysisDriver
from analysis_driver.pipelines.variant_calling import GATKStage, BaseRecal, PrintReads, \
    RealignTarget, Realign, HaplotypeCaller, GenotypeGVCFs, SelectVariants, VariantFiltration
from unittest.mock import patch, call

patch_executor = patch('analysis_driver.pipelines.variant_calling.executor.execute')


class TestVariantCalling(TestAnalysisDriver):
    stage_cls = GATKStage
    extra_args = {}
    dataset = NamedMock(
        real_name='test_dataset',
        reference_genome='reference_genome',
        user_sample_id='test_user_sample_id',
        genome_version='genome_version'
    )

    def setUp(self):
        self.stage = self.stage_cls(dataset=self.dataset, **self.extra_args)


class TestGATKGenericStage(TestVariantCalling):
    def test_gatk_run_dir(self):
        run_dir = self.stage.gatk_run_dir
        assert run_dir == 'tests/assets/jobs/test_dataset/gatk_var_calling'
        assert os.path.isdir('tests/assets/jobs/test_dataset/gatk_var_calling')

    def test_gatk_cmd(self):
        minimal_cmd = self.stage.gatk_cmd('GATKClass', 'test_outfile')
        assert minimal_cmd == ('java -Djava.io.tmpdir=tests/assets/jobs/test_dataset/gatk_var_calling -XX:+UseSerialGC '
                               '-Xmx16G -jar path/to/gatk -R reference_genome -T GATKClass --read_filter BadCigar '
                               '--read_filter NotPrimaryAlignment -o test_outfile -l INFO -U LENIENT_VCF_PROCESSING '
                               '-nct 16 -nt 16')

        full_cmd = self.stage.gatk_cmd('GATKClass', 'test_outfile', 'test_bam', 17, 18, 19, ' --flag option')
        assert full_cmd == ('java -Djava.io.tmpdir=tests/assets/jobs/test_dataset/gatk_var_calling -XX:+UseSerialGC '
                            '-Xmx17G -jar path/to/gatk -R reference_genome -T GATKClass --read_filter BadCigar '
                            '--read_filter NotPrimaryAlignment -o test_outfile -l INFO -U LENIENT_VCF_PROCESSING '
                            '-I test_bam --flag option -nct 18 -nt 19')

    def test_basename(self):
        basename = self.stage.basename
        assert basename == 'tests/assets/jobs/test_dataset/gatk_var_calling/test_user_sample_id'
        os.path.isdir(basename)

    def test_sorted_bam(self):
        assert self.stage.sorted_bam == 'tests/assets/jobs/test_dataset/test_dataset.bam'

    def test_recal_bam(self):
        assert self.stage.recal_bam == 'tests/assets/jobs/test_dataset/gatk_var_calling/test_user_sample_id_recal.bam'

    def test_output_grp(self):
        assert self.stage.output_grp == 'tests/assets/jobs/test_dataset/gatk_var_calling/test_user_sample_id.grp'

    def test_output_intervals(self):
        assert self.stage.output_intervals == 'tests/assets/jobs/test_dataset/gatk_var_calling/test_user_sample_id.intervals'

    def test_indel_realigned_bam(self):
        assert self.stage.indel_realigned_bam == 'tests/assets/jobs/test_dataset/gatk_var_calling/test_user_sample_id_indel_realigned.bam'

    def test_dataset_gvcf(self):
        assert self.stage.sample_gvcf == 'tests/assets/jobs/test_dataset/gatk_var_calling/test_user_sample_id.g.vcf'

    def test_genotyped_vcf(self):
        assert self.stage.genotyped_vcf == 'tests/assets/jobs/test_dataset/gatk_var_calling/test_user_sample_id.vcf'

    def test_raw_snp_vcf(self):
        assert self.stage.raw_snp_vcf == 'tests/assets/jobs/test_dataset/gatk_var_calling/test_user_sample_id_raw_snp.vcf'

    def test_filter_snp_vcf(self):
        assert self.stage.filter_snp_vcf == 'tests/assets/jobs/test_dataset/gatk_var_calling/test_user_sample_id_filter_snp.vcf'

    def test_dbsnp(self):
        assert self.stage.dbsnp == '/path/to/dbsnp.vcf.gz'

    def test_known_indels(self):
        assert self.stage.known_indels == '/path/to/known/indels'


class TestGATKStage(TestVariantCalling):
    def setUp(self):
        super().setUp()
        self.patched_executor = patch('analysis_driver.pipelines.variant_calling.executor.execute')
        self.mocked_executor = self.patched_executor.start()
        self.patched_gatk_cmd = patch.object(self.stage_cls, 'gatk_cmd', return_value='a_cmd')
        self.mocked_gatk_cmd = self.patched_gatk_cmd.start()

    def tearDown(self):
        self.patched_executor.stop()
        self.patched_gatk_cmd.stop()

    def _test_bgzip_and_tabix(self, vcf_file):
        assert self.mocked_executor.call_args_list[1] == call(
            'path/to/bgzip -f ' + vcf_file, cpus=1, job_name='bgzip', mem=8,
            working_dir='tests/assets/jobs/test_dataset/gatk_var_calling'
        )
        assert self.mocked_executor.call_args_list[2] == call(
            'path/to/tabix -f -p vcf ' + vcf_file + '.gz', cpus=1, job_name='tabix', mem=8,
            working_dir='tests/assets/jobs/test_dataset/gatk_var_calling'
        )


class TestBaseRecal(TestGATKStage):
    stage_cls = BaseRecal

    def test_run(self):
        self.stage._run()
        assert self.mocked_executor.call_count == 1
        self.mocked_gatk_cmd.assert_called_with(
            'BaseRecalibrator',
            'tests/assets/jobs/test_dataset/gatk_var_calling/test_user_sample_id.grp',
            input_bam='tests/assets/jobs/test_dataset/test_dataset.bam',
            xmx=48,
            nt=1,
            ext=' --knownSites /path/to/dbsnp.vcf.gz'
        )
        self.mocked_executor.assert_called_with(
            'a_cmd',
            cpus=16,
            job_name='gatk_base_recal',
            mem=64,
            working_dir='tests/assets/jobs/test_dataset/gatk_var_calling'
        )


class TestPrintReads(TestGATKStage):
    stage_cls = PrintReads

    def test_run(self):
        self.stage._run()
        assert self.mocked_executor.call_count == 1
        self.mocked_gatk_cmd.assert_called_with(
            'PrintReads',
            'tests/assets/jobs/test_dataset/gatk_var_calling/test_user_sample_id_recal.bam',
            input_bam='tests/assets/jobs/test_dataset/test_dataset.bam',
            xmx=48,
            nt=1,
            ext=' -BQSR tests/assets/jobs/test_dataset/gatk_var_calling/test_user_sample_id.grp'
        )
        self.mocked_executor.assert_called_with(
            'a_cmd',
            cpus=16,
            job_name='gatk_print_reads',
            mem=64,
            working_dir='tests/assets/jobs/test_dataset/gatk_var_calling'
        )


class TestRealignTarget(TestGATKStage):
    stage_cls = RealignTarget

    def test_run(self):
        self.stage._run()
        assert self.mocked_executor.call_count == 1
        self.mocked_gatk_cmd.assert_called_with(
            'RealignerTargetCreator',
            'tests/assets/jobs/test_dataset/gatk_var_calling/test_user_sample_id.intervals',
            input_bam='tests/assets/jobs/test_dataset/gatk_var_calling/test_user_sample_id_recal.bam',
            xmx=32,
            nct=1,
            nt=1,
            ext=' --known /path/to/known/indels'
        )
        self.mocked_executor.assert_called_with(
            'a_cmd',
            job_name='gatk_realign_target',
            mem=32,
            working_dir='tests/assets/jobs/test_dataset/gatk_var_calling'
        )


class TestRealign(TestGATKStage):
    stage_cls = Realign

    def test_run(self):
        self.stage._run()
        assert self.mocked_executor.call_count == 1
        self.mocked_gatk_cmd.assert_called_with(
            'IndelRealigner',
            'tests/assets/jobs/test_dataset/gatk_var_calling/test_user_sample_id_indel_realigned.bam',
            input_bam='tests/assets/jobs/test_dataset/gatk_var_calling/test_user_sample_id_recal.bam',
            xmx=32,
            nct=1,
            nt=1,
            ext=' -targetIntervals tests/assets/jobs/test_dataset/gatk_var_calling/test_user_sample_id.intervals '
                '--knownAlleles /path/to/known/indels'
        )
        self.mocked_executor.assert_called_with(
            'a_cmd',
            job_name='gatk_indel_realign',
            mem=32,
            working_dir='tests/assets/jobs/test_dataset/gatk_var_calling'
        )


class TestHaplotypeCaller(TestGATKStage):
    stage_cls = HaplotypeCaller
    extra_args = {'input_bam': 'test_bam'}

    def test_run(self):
        self.stage._run()
        assert self.mocked_executor.call_count == 3  # Command + bgzip + tabix
        self.mocked_gatk_cmd.assert_called_with(
            'HaplotypeCaller',
            'tests/assets/jobs/test_dataset/gatk_var_calling/test_user_sample_id.g.vcf',
            input_bam='test_bam',
            xmx=48,
            nt=1,
            ext=' --pair_hmm_implementation VECTOR_LOGLESS_CACHING '
                '-ploidy 2 '
                '--emitRefConfidence GVCF '
                '--variant_index_type LINEAR '
                '--variant_index_parameter 128000  '
                '--annotation BaseQualityRankSumTest '
                '--annotation FisherStrand '
                '--annotation GCContent '
                '--annotation HaplotypeScore '
                '--annotation HomopolymerRun '
                '--annotation MappingQualityRankSumTest '
                '--annotation MappingQualityZero '
                '--annotation QualByDepth '
                '--annotation ReadPosRankSumTest '
                '--annotation RMSMappingQuality '
                '--annotation DepthPerAlleleBySample '
                '--annotation Coverage '
                '--annotation ClippingRankSumTest '
                '--annotation DepthPerSampleHC '
                '--dbsnp /path/to/dbsnp.vcf.gz'
        )
        assert self.mocked_executor.call_args_list[0] == call(
            'a_cmd',
            cpus=16,
            job_name='gatk_haplotype_call',
            mem=64,
            working_dir='tests/assets/jobs/test_dataset/gatk_var_calling'
        )
        self._test_bgzip_and_tabix('tests/assets/jobs/test_dataset/gatk_var_calling/test_user_sample_id.g.vcf')


class TestGenotypeGVCFs(TestGATKStage):
    stage_cls = GenotypeGVCFs

    def test_run(self):
        self.stage._run()

        assert self.mocked_executor.call_count == 3  # Command + bgzip + tabix
        self.mocked_gatk_cmd.assert_called_with(
            'GenotypeGVCFs',
            'tests/assets/jobs/test_dataset/gatk_var_calling/test_user_sample_id.vcf',
            nct=1,
            ext=' --variant tests/assets/jobs/test_dataset/gatk_var_calling/test_user_sample_id.g.vcf.gz'
        )
        assert self.mocked_executor.call_args_list[0] == call(
            'a_cmd',
            job_name='gatk_genotype_gvcfs',
            mem=16,
            working_dir='tests/assets/jobs/test_dataset/gatk_var_calling'
        )
        self._test_bgzip_and_tabix('tests/assets/jobs/test_dataset/gatk_var_calling/test_user_sample_id.vcf')


class TestSelectVariants(TestGATKStage):
    stage_cls = SelectVariants

    def test_run(self):
        self.stage._run()

        assert self.mocked_executor.call_count == 3  # Command + bgzip + tabix
        self.mocked_gatk_cmd.assert_called_with(
            'SelectVariants',
            'tests/assets/jobs/test_dataset/gatk_var_calling/test_user_sample_id_raw_snp.vcf',
            nct=1,
            nt=16,
            ext=' -V tests/assets/jobs/test_dataset/gatk_var_calling/test_user_sample_id.vcf.gz -selectType SNP'
        )
        assert self.mocked_executor.call_args_list[0] == call(
            'a_cmd',
            job_name='var_filtration',
            mem=16,
            working_dir='tests/assets/jobs/test_dataset/gatk_var_calling'
        )
        self._test_bgzip_and_tabix('tests/assets/jobs/test_dataset/gatk_var_calling/test_user_sample_id_raw_snp.vcf')


class TestVariantFiltration(TestGATKStage):
    stage_cls = VariantFiltration

    def test_run(self):
        self.stage._run()
        assert self.mocked_executor.call_count == 3  # Command + bgzip + tabix
        self.mocked_gatk_cmd.assert_called_with(
            'VariantFiltration',
            'tests/assets/jobs/test_dataset/gatk_var_calling/test_user_sample_id_filter_snp.vcf',
            nct=1,
            nt=1,
            ext=' -V tests/assets/jobs/test_dataset/gatk_var_calling/test_user_sample_id_raw_snp.vcf.gz '
                "--filterExpression 'QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0' "
                "--filterName 'SNP_FILTER'"
        )
        assert self.mocked_executor.call_args_list[0] == call(
            'a_cmd',
            job_name='var_filtration',
            mem=16,
            working_dir='tests/assets/jobs/test_dataset/gatk_var_calling'
        )
        self._test_bgzip_and_tabix('tests/assets/jobs/test_dataset/gatk_var_calling/test_user_sample_id_filter_snp.vcf')
