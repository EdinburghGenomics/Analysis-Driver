import os
from unittest.mock import patch
from tests.test_analysisdriver import NamedMock, TestAnalysisDriver
from analysis_driver.pipelines import variant_calling


fake_dataset = NamedMock(
    real_name='test_sample',
    reference_genome='reference_genome',
    user_sample_id='test_user_sample_id',
    genome_version='genome_version'
)


class TestVarCallingStage(TestAnalysisDriver):
    cls = variant_calling.GATKStage
    exp_command = None
    exec_kwargs = {}
    exp_run_dir = 'tests/assets/jobs/test_sample/gatk_var_calling'

    def setUp(self):
        self.g = self.cls(dataset=fake_dataset)
        self.patched_java_cmd = patch('analysis_driver.pipelines.variant_calling.java_command', return_value='[java] ')
        self.mocked_java_cmd = self.patched_java_cmd.start()

    def tearDown(self):
        self.mocked_java_cmd = self.patched_java_cmd.stop()

    def test_gatk_run_dir(self):
        assert self.g.gatk_run_dir == self.exp_run_dir
        assert os.path.isdir(self.exp_run_dir)

    def test_gatk_cmd(self):
        minimal_gatk_cmd = self.g.gatk_cmd('GATKfunction', 'test_outfile')
        assert minimal_gatk_cmd == (
            '[java] -R reference_genome -T GATKfunction --read_filter BadCigar --read_filter NotPrimaryAlignment -o '
            'test_outfile -l INFO -U LENIENT_VCF_PROCESSING -nct 16 -nt 16'
        )
        self.mocked_java_cmd.assert_called_with(memory=16, tmp_dir=self.exp_run_dir, jar='path/to/gatk')
        self.mocked_java_cmd.reset_mock()

        gatk_cmd = self.g.gatk_cmd('GATKfunction', 'test_outfile', 'test_bam', ext=' --flag option')
        assert gatk_cmd == (
            '[java] -R reference_genome -T GATKfunction --read_filter BadCigar --read_filter NotPrimaryAlignment -o '
            'test_outfile -l INFO -U LENIENT_VCF_PROCESSING -I test_bam --flag option -nct 16 -nt 16'
        )
        self.mocked_java_cmd.assert_called_with(memory=16, tmp_dir=self.exp_run_dir, jar='path/to/gatk')

    def test_basename(self):
        assert self.g.basename == os.path.join(self.exp_run_dir, 'test_user_sample_id')

    def test_sorted_bam(self):
        assert self.g.sorted_bam == 'tests/assets/jobs/test_sample/test_sample.bam'

    def test_recal_bam(self):
        assert self.g.recal_bam == os.path.join(self.exp_run_dir, 'test_user_sample_id_recal.bam')

    def test_output_grp(self):
        assert self.g.output_grp == os.path.join(self.exp_run_dir, 'test_user_sample_id.grp')

    def test_output_intervals(self):
        assert self.g.output_intervals == os.path.join(self.exp_run_dir, 'test_user_sample_id.intervals')

    def test_indel_realigned_bam(self):
        assert self.g.indel_realigned_bam == os.path.join(self.exp_run_dir, 'test_user_sample_id_indel_realigned.bam')

    def test_sample_gvcf(self):
        assert self.g.sample_gvcf == os.path.join(self.exp_run_dir, 'test_user_sample_id.g.vcf')

    def test_genotyped_vcf(self):
        assert self.g.genotyped_vcf == os.path.join(self.exp_run_dir, 'test_user_sample_id.vcf')

    def test_raw_snp_vcf(self):
        assert self.g.raw_snp_vcf == os.path.join(self.exp_run_dir, 'test_user_sample_id_raw_snp.vcf')

    def test_filter_snp_vcf(self):
        assert self.g.filter_snp_vcf == os.path.join(self.exp_run_dir, 'test_user_sample_id_filter_snp.vcf')

    def test_dbsnp(self):
        assert self.g.dbsnp == '/path/to/dbsnp.vcf.gz'

    def test_known_indels(self):
        assert self.g.known_indels == '/path/to/known/indels'

    def _check_calls(self, mocked_execute):
        assert mocked_execute.call_count == 1
        positional_args = self.exp_command
        if isinstance(positional_args, str):
            positional_args = (positional_args,)

        mocked_execute.assert_called_with(
            *positional_args,
            working_dir=self.exp_run_dir,
            **self.exec_kwargs
        )

    @patch('analysis_driver.pipelines.variant_calling.executor.execute')
    def test_run(self, mocked_execute):
        return_value = self.g._run()
        if return_value:  # mocked return value from patched execute
            self._check_calls(mocked_execute)


class TestGATKStage(TestVarCallingStage):
    xmx = None

    def _check_calls(self, mocked_execute):
        super()._check_calls(mocked_execute)
        self.mocked_java_cmd.assert_called_with(
            memory=self.xmx,
            tmp_dir=self.exp_run_dir,
            jar='path/to/gatk'
        )


class TestBaseRecal(TestGATKStage):
    cls = variant_calling.BaseRecal
    exp_command = (
        '[java] '
        '-R reference_genome '
        '-T BaseRecalibrator '
        '--read_filter BadCigar '
        '--read_filter NotPrimaryAlignment '
        '-o tests/assets/jobs/test_sample/gatk_var_calling/test_user_sample_id.grp '
        '-l INFO '
        '-U LENIENT_VCF_PROCESSING '
        '-I tests/assets/jobs/test_sample/test_sample.bam '
        '--knownSites /path/to/dbsnp.vcf.gz '
        '-nct 16'
    )
    xmx = 48
    exec_kwargs = {'cpus': 16, 'job_name': 'gatk_base_recal', 'mem': 64}


class TestPrintReads(TestGATKStage):
    cls = variant_calling.PrintReads
    exp_command = (
        '[java] '
        '-R reference_genome '
        '-T PrintReads '
        '--read_filter BadCigar '
        '--read_filter NotPrimaryAlignment '
        '-o tests/assets/jobs/test_sample/gatk_var_calling/test_user_sample_id_recal.bam '
        '-l INFO -U LENIENT_VCF_PROCESSING '
        '-I tests/assets/jobs/test_sample/test_sample.bam '
        '-BQSR tests/assets/jobs/test_sample/gatk_var_calling/test_user_sample_id.grp '
        '-nct 16'
    )
    xmx = 48
    exec_kwargs = {'cpus': 16, 'job_name': 'gatk_print_reads', 'mem': 64}


class TestRealignTarget(TestGATKStage):
    cls = variant_calling.RealignTarget
    exp_command = (
        '[java] '
        '-R reference_genome '
        '-T RealignerTargetCreator '
        '--read_filter BadCigar '
        '--read_filter NotPrimaryAlignment '
        '-o tests/assets/jobs/test_sample/gatk_var_calling/test_user_sample_id.intervals '
        '-l INFO '
        '-U LENIENT_VCF_PROCESSING '
        '-I tests/assets/jobs/test_sample/gatk_var_calling/test_user_sample_id_recal.bam '
        '--known /path/to/known/indels'
    )
    xmx = 16
    exec_kwargs = {'job_name': 'gatk_realign_target', 'mem': 16}


class TestRealign(TestGATKStage):
    cls = variant_calling.Realign
    exp_command = (
        '[java] '
        '-R reference_genome '
        '-T IndelRealigner '
        '--read_filter BadCigar '
        '--read_filter NotPrimaryAlignment '
        '-o tests/assets/jobs/test_sample/gatk_var_calling/test_user_sample_id_indel_realigned.bam '
        '-l INFO -U LENIENT_VCF_PROCESSING '
        '-I tests/assets/jobs/test_sample/gatk_var_calling/test_user_sample_id_recal.bam '
        '-targetIntervals tests/assets/jobs/test_sample/gatk_var_calling/test_user_sample_id.intervals '
        '--knownAlleles /path/to/known/indels'
    )
    xmx = 16
    exec_kwargs = {'job_name': 'gatk_indel_realign', 'mem': 16}


class TestHaplotypeCaller(TestGATKStage):
    exp_command = (
        '[java] '
        '-R reference_genome '
        '-T HaplotypeCaller '
        '--read_filter BadCigar '
        '--read_filter NotPrimaryAlignment '
        '-o tests/assets/jobs/test_sample/gatk_var_calling/test_user_sample_id.g.vcf '
        '-l INFO '
        '-U LENIENT_VCF_PROCESSING '
        '-I test_bam --pair_hmm_implementation VECTOR_LOGLESS_CACHING '
        '-ploidy 2 '
        '--emitRefConfidence GVCF '
        '--variant_index_type LINEAR '
        '--variant_index_parameter 128000  '
        '-nct 16 '
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
    xmx = 48
    exec_kwargs = {'cpus': 16, 'job_name': 'gatk_haplotype_call', 'mem': 64}

    def setUp(self):
        super().setUp()
        self.g = variant_calling.HaplotypeCaller(dataset=fake_dataset, input_bam='test_bam')


class TestGenotypeGVCFs(TestGATKStage):
    cls = variant_calling.GenotypeGVCFs
    exp_command = (
        '[java] '
        '-R reference_genome '
        '-T GenotypeGVCFs '
        '--read_filter BadCigar '
        '--read_filter NotPrimaryAlignment '
        '-o tests/assets/jobs/test_sample/gatk_var_calling/test_user_sample_id.vcf '
        '-l INFO '
        '-U LENIENT_VCF_PROCESSING '
        '--variant tests/assets/jobs/test_sample/gatk_var_calling/test_user_sample_id.g.vcf '
        '-nt 16'
    )
    xmx = 16
    exec_kwargs = {'job_name': 'gatk_genotype_gvcfs', 'mem': 16}


class TestSelectVariants(TestGATKStage):
    cls = variant_calling.SelectVariants
    exp_command = (
        '[java] '
        '-R reference_genome '
        '-T SelectVariants '
        '--read_filter BadCigar '
        '--read_filter NotPrimaryAlignment '
        '-o tests/assets/jobs/test_sample/gatk_var_calling/test_user_sample_id_raw_snp.vcf '
        '-l INFO -U LENIENT_VCF_PROCESSING '
        '-nt 16 '
        '-V tests/assets/jobs/test_sample/gatk_var_calling/test_user_sample_id.vcf '
        '-selectType SNP '
    )
    xmx = 16
    exec_kwargs = {'job_name': 'var_filtration', 'mem': 16}


class TestVariantFiltration(TestGATKStage):
    cls = variant_calling.VariantFiltration
    exp_command = (
        '[java] '
        '-R reference_genome '
        '-T VariantFiltration '
        '--read_filter BadCigar '
        '--read_filter NotPrimaryAlignment '
        '-o tests/assets/jobs/test_sample/gatk_var_calling/test_user_sample_id_filter_snp.vcf '
        '-l INFO '
        '-U LENIENT_VCF_PROCESSING '
        '-V tests/assets/jobs/test_sample/gatk_var_calling/test_user_sample_id_raw_snp.vcf '
        "--filterExpression 'QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0' "
        "--filterName 'SNP_FILTER'"
    )
    xmx = 16
    exec_kwargs = {'job_name': 'var_filtration', 'mem': 16}


class TestBGZip(TestVarCallingStage):
    cls = variant_calling.BGZip
    exp_command = (
        'path/to/bgzip ' + os.path.join(TestVarCallingStage.exp_run_dir, 'test_user_sample_id.g.vcf'),
        'path/to/bgzip ' + os.path.join(TestVarCallingStage.exp_run_dir, 'test_user_sample_id_filter_snp.vcf')
    )
    exec_kwargs = {'cpus': 1, 'job_name': 'bgzip', 'mem': 8}


class TestTabix(TestVarCallingStage):
    cls = variant_calling.Tabix
    exp_command = (
        'path/to/tabix -p vcf ' + os.path.join(TestVarCallingStage.exp_run_dir, 'test_user_sample_id.g.vcf.gz'),
        'path/to/tabix -p vcf ' + os.path.join(TestVarCallingStage.exp_run_dir, 'test_user_sample_id_filter_snp.vcf.gz'),
    )
    exec_kwargs = {'cpus': 1, 'job_name': 'tabix', 'mem': 8}
