import os
from tests.test_analysisdriver import NamedMock, TestAnalysisDriver
from analysis_driver.pipelines.variant_calling import GATKStage, BaseRecal, PrintReads, \
    RealignTarget, Realign, HaplotypeCaller, GenotypeGVCFs, SelectVariants, VariantFiltration
from unittest.mock import patch, call

patch_executor = patch('analysis_driver.pipelines.variant_calling.executor.execute')


class TestGATKStage():
    dataset = NamedMock(real_name='test_sample',
                        reference_genome='test_reference',
                        user_sample_id='test_user_sample_id',
                        genome_version='genome_version',
                        genome_dict={
                            'data_files': {
                                'variation': 'path/to/genomes_dir/path/to/dbsnp.vcf.gz'
                            }
                        })

    g = GATKStage(dataset=dataset)

    def test_gatk_run_dir(self):
        run_dir = self.g.gatk_run_dir
        assert run_dir == 'tests/assets/jobs/test_sample/gatk_var_calling'
        assert os.path.isdir('tests/assets/jobs/test_sample/gatk_var_calling')

    def test_gatk_cmd(self):
        gatk_cmd_without_bam = self.g.gatk_cmd('GATKfunction', 'test_outfile')
        assert gatk_cmd_without_bam == 'path/to/java_8 -Djava.io.tmpdir=tests/assets/jobs/test_sample/gatk_var_calling ' \
                           '-XX:+UseSerialGC -Xmx16G -jar path/to/gatk_4 -R test_reference -T GATKfunction ' \
                           '--read_filter BadCigar --read_filter NotPrimaryAlignment -o test_outfile ' \
                           '-l INFO -U LENIENT_VCF_PROCESSING -nct 16 -nt 16'

        gatk_cmd_with_bam = self.g.gatk_cmd('GATKfunction', 'test_outfile', 'test_bam')
        assert gatk_cmd_with_bam == 'path/to/java_8 -Djava.io.tmpdir=tests/assets/jobs/test_sample/gatk_var_calling ' \
                           '-XX:+UseSerialGC -Xmx16G -jar path/to/gatk_4 -R test_reference -T GATKfunction ' \
                           '--read_filter BadCigar --read_filter NotPrimaryAlignment -o test_outfile ' \
                           '-l INFO -U LENIENT_VCF_PROCESSING -I test_bam -nct 16 -nt 16'

        gatk_option = 'option'
        gatk_cmd_with_ext = self.g.gatk_cmd('GATKfunction', 'test_outfile', ext=' --flag ' + gatk_option)
        assert gatk_cmd_with_ext == 'path/to/java_8 -Djava.io.tmpdir=tests/assets/jobs/test_sample/gatk_var_calling ' \
                           '-XX:+UseSerialGC -Xmx16G -jar path/to/gatk_4 -R test_reference -T GATKfunction ' \
                           '--read_filter BadCigar --read_filter NotPrimaryAlignment -o test_outfile ' \
                           '-l INFO -U LENIENT_VCF_PROCESSING --flag option -nct 16 -nt 16'

    def test_basename(self):
        basename = self.g.basename
        assert basename == 'tests/assets/jobs/test_sample/gatk_var_calling/test_user_sample_id'
        os.path.isdir(basename)

    def test_sorted_bam(self):
        bam = self.g.sorted_bam
        assert bam == 'tests/assets/jobs/test_sample/test_sample.bam'

    def test_recal_bam(self):
        recal_bam = self.g.recal_bam
        assert recal_bam == 'tests/assets/jobs/test_sample/gatk_var_calling/test_user_sample_id_recal.bam'

    def test_output_grp(self):
        output_grp = self.g.output_grp
        assert output_grp == 'tests/assets/jobs/test_sample/gatk_var_calling/test_user_sample_id.grp'

    def test_output_intervals(self):
        intervals = self.g.output_intervals
        assert intervals == 'tests/assets/jobs/test_sample/gatk_var_calling/test_user_sample_id.intervals'

    def test_indel_realigned_bam(self):
        indel_realigned_bam = self.g.indel_realigned_bam
        assert indel_realigned_bam == 'tests/assets/jobs/test_sample/gatk_var_calling/test_user_sample_id_indel_realigned.bam'

    def test_sample_gvcf(self):
        gvcf = self.g.sample_gvcf
        assert gvcf == 'tests/assets/jobs/test_sample/gatk_var_calling/test_user_sample_id.g.vcf'

    def test_genotyped_vcf(self):
        genotyped_vcf = self.g.genotyped_vcf
        assert genotyped_vcf == 'tests/assets/jobs/test_sample/gatk_var_calling/test_user_sample_id.vcf'

    def test_raw_snp_vcf(self):
        raw_snp_vcf = self.g.raw_snp_vcf
        assert raw_snp_vcf == 'tests/assets/jobs/test_sample/gatk_var_calling/test_user_sample_id_raw_snp.vcf'

    def test_filter_snp_vcf(self):
        filter_snp_vcf = self.g.filter_snp_vcf
        assert filter_snp_vcf == 'tests/assets/jobs/test_sample/gatk_var_calling/test_user_sample_id_filter_snp.vcf'

    def test_dbsnp(self):
        dbsnp = self.g.dbsnp
        assert dbsnp == 'path/to/genomes_dir/path/to/dbsnp.vcf.gz'


class TestVariantCalling(TestAnalysisDriver):

    def setUp(self):
        self.dataset = NamedMock(
            real_name='test_dataset',
            user_sample_id='test_user_sample_id',
            genome_version='genome_version',
            reference_genome='reference_genome',
            genome_dict={
                'data_files': {
                    'variation': 'path/to/genomes_dir/path/to/dbsnp.vcf.gz',
                    'vqsr': {
                        'hapmap': 'path/to/hapmap_annotation',
                        'omni': 'path/to/omni_annotation',
                        'thousand_genomes': 'path/to/1000g_annotation',
                        'dbsnp': 'path/to/dbsnp_annotation',
                        'mills': 'path/to/mills_annotation'
                    },
                }
            })


class TestBaseRecal(TestVariantCalling):

    def test_run(self):
        b = BaseRecal(dataset=self.dataset)
        with patch_executor as e:
            b._run()
            assert e.call_count == 1
            e.assert_called_with("path/to/java_8 -Djava.io.tmpdir=tests/assets/jobs/test_dataset/gatk_var_calling "
                                 "-XX:+UseSerialGC "
                                 "-Xmx48G "
                                 "-jar path/to/gatk_4 "
                                 "-R reference_genome "
                                 "-T BaseRecalibrator "
                                 "--read_filter BadCigar "
                                 "--read_filter NotPrimaryAlignment "
                                 "-o tests/assets/jobs/test_dataset/gatk_var_calling/test_user_sample_id.grp "
                                 "-l INFO "
                                 "-U LENIENT_VCF_PROCESSING "
                                 "-I tests/assets/jobs/test_dataset/test_dataset.bam "
                                 "--knownSites path/to/genomes_dir/path/to/dbsnp.vcf.gz "
                                 "-nct 16",
                                 cpus=16,
                                 job_name='gatk_base_recal',
                                 mem=64,
                                 working_dir='tests/assets/jobs/test_dataset/gatk_var_calling')


class TestPrintReads(TestVariantCalling):

    def test_run(self):
        p = PrintReads(dataset=self.dataset)

        with patch_executor as e:
            p._run()
            assert e.call_count == 1
            e.assert_called_with('path/to/java_8 -Djava.io.tmpdir=tests/assets/jobs/test_dataset/gatk_var_calling '
                                 '-XX:+UseSerialGC '
                                 '-Xmx48G '
                                 '-jar path/to/gatk_4 '
                                 '-R reference_genome '
                                 '-T PrintReads '
                                 '--read_filter BadCigar '
                                 '--read_filter NotPrimaryAlignment '
                                 '-o tests/assets/jobs/test_dataset/gatk_var_calling/test_user_sample_id_recal.bam '
                                 '-l INFO -U LENIENT_VCF_PROCESSING '
                                 '-I tests/assets/jobs/test_dataset/test_dataset.bam '
                                 '-BQSR tests/assets/jobs/test_dataset/gatk_var_calling/test_user_sample_id.grp '
                                 '-nct 16',
                                 cpus=16,
                                 job_name='gatk_print_reads',
                                 mem=64,
                                 working_dir='tests/assets/jobs/test_dataset/gatk_var_calling')


class TestRealignTarget(TestVariantCalling):

    def test_run(self):
        p = RealignTarget(dataset=self.dataset)

        with patch_executor as e:
            p._run()
            assert e.call_count == 1
            e.assert_called_with('path/to/java_8 -Djava.io.tmpdir=tests/assets/jobs/test_dataset/gatk_var_calling '
                                 '-XX:+UseSerialGC '
                                 '-Xmx32G '
                                 '-jar path/to/gatk_4 '
                                 '-R reference_genome '
                                 '-T RealignerTargetCreator '
                                 '--read_filter BadCigar '
                                 '--read_filter NotPrimaryAlignment '
                                 '-o tests/assets/jobs/test_dataset/gatk_var_calling/test_user_sample_id.intervals '
                                 '-l INFO '
                                 '-U LENIENT_VCF_PROCESSING '
                                 '-I tests/assets/jobs/test_dataset/gatk_var_calling/test_user_sample_id_recal.bam',
                                 job_name='gatk_realign_target',
                                 mem=32,
                                 working_dir='tests/assets/jobs/test_dataset/gatk_var_calling')


class TestRealign(TestVariantCalling):

    def test_run(self):
        p = Realign(dataset=self.dataset)
        with patch_executor as e:
            p._run()
            assert e.call_count == 1
            e.assert_called_with('path/to/java_8 -Djava.io.tmpdir=tests/assets/jobs/test_dataset/gatk_var_calling '
                                 '-XX:+UseSerialGC '
                                 '-Xmx32G '
                                 '-jar path/to/gatk_4 '
                                 '-R reference_genome '
                                 '-T IndelRealigner'
                                 ' --read_filter BadCigar '
                                 '--read_filter NotPrimaryAlignment '
                                 '-o tests/assets/jobs/test_dataset/gatk_var_calling/test_user_sample_id_indel_realigned.bam '
                                 '-l INFO -U LENIENT_VCF_PROCESSING '
                                 '-I tests/assets/jobs/test_dataset/gatk_var_calling/test_user_sample_id_recal.bam '
                                 '-targetIntervals tests/assets/jobs/test_dataset/gatk_var_calling/test_user_sample_id.intervals',
                                 job_name='gatk_indel_realign',
                                 mem=32,
                                 working_dir='tests/assets/jobs/test_dataset/gatk_var_calling')


class TestHaplotypeCaller(TestVariantCalling):

    def test_run(self):
        p = HaplotypeCaller(dataset=self.dataset, input_bam='test_bam')

        with patch_executor as e:
            p._run()
            assert e.call_count == 3  # Command + bgzip + tabix
            assert e.call_args_list[0] == call(
                'path/to/java_8 -Djava.io.tmpdir=tests/assets/jobs/test_dataset/gatk_var_calling '
                '-XX:+UseSerialGC '
                '-Xmx48G '
                '-jar path/to/gatk_4 '
                '-R reference_genome '
                '-T HaplotypeCaller '
                '--read_filter BadCigar '
                '--read_filter NotPrimaryAlignment '
                '-o tests/assets/jobs/test_dataset/gatk_var_calling/test_user_sample_id.g.vcf '
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
                '--dbsnp path/to/genomes_dir/path/to/dbsnp.vcf.gz',
                cpus=16,
                job_name='gatk_haplotype_call',
                mem=64,
                working_dir='tests/assets/jobs/test_dataset/gatk_var_calling'
            )
            _test_bgzip_and_tabix(e, 'tests/assets/jobs/test_dataset/gatk_var_calling/test_user_sample_id.g.vcf')


class TestGenotypeGVCFs(TestVariantCalling):

    def test_run(self):
        p = GenotypeGVCFs(dataset=self.dataset)

        with patch_executor as e:
            p._run()
            assert e.call_count == 3  # Command + bgzip + tabix
            assert e.call_args_list[0] == call(
                'path/to/java_8 -Djava.io.tmpdir=tests/assets/jobs/test_dataset/gatk_var_calling '
                '-XX:+UseSerialGC '
                '-Xmx16G '
                '-jar path/to/gatk_4 '
                '-R reference_genome '
                '-T GenotypeGVCFs '
                '--read_filter BadCigar '
                '--read_filter NotPrimaryAlignment '
                '-o tests/assets/jobs/test_dataset/gatk_var_calling/test_user_sample_id.vcf '
                '-l INFO '
                '-U LENIENT_VCF_PROCESSING '
                '--variant tests/assets/jobs/test_dataset/gatk_var_calling/test_user_sample_id.g.vcf.gz '
                '-nt 16',
                job_name='gatk_genotype_gvcfs',
                mem=16,
                working_dir='tests/assets/jobs/test_dataset/gatk_var_calling'
            )
            _test_bgzip_and_tabix(e, 'tests/assets/jobs/test_dataset/gatk_var_calling/test_user_sample_id.vcf')


class TestSelectVariants(TestVariantCalling):

    def test_run(self):
        p = SelectVariants(dataset=self.dataset)

        with patch_executor as e:
            p._run()
            assert e.call_count == 3  # Command + bgzip + tabix
            assert e.call_args_list[0] == call(
                'path/to/java_8 -Djava.io.tmpdir=tests/assets/jobs/test_dataset/gatk_var_calling '
                '-XX:+UseSerialGC '
                '-Xmx16G '
                '-jar path/to/gatk_4 '
                '-R reference_genome '
                '-T SelectVariants '
                '--read_filter BadCigar '
                '--read_filter NotPrimaryAlignment '
                '-o tests/assets/jobs/test_dataset/gatk_var_calling/test_user_sample_id_raw_snp.vcf '
                '-l INFO -U LENIENT_VCF_PROCESSING '
                '-nt 16 '
                '-V tests/assets/jobs/test_dataset/gatk_var_calling/test_user_sample_id.vcf.gz '
                '-selectType SNP ',
                job_name='var_filtration',
                mem=16,
                working_dir='tests/assets/jobs/test_dataset/gatk_var_calling'
            )
            _test_bgzip_and_tabix(e, 'tests/assets/jobs/test_dataset/gatk_var_calling/test_user_sample_id_raw_snp.vcf')


class TestVariantFiltration(TestVariantCalling):

    def test_run(self):
        p = VariantFiltration(dataset=self.dataset)

        with patch_executor as e:
            p._run()
            assert e.call_count == 3  # Command + bgzip + tabix
            assert e.call_args_list[0] == call(
                "path/to/java_8 -Djava.io.tmpdir=tests/assets/jobs/test_dataset/gatk_var_calling "
                "-XX:+UseSerialGC "
                "-Xmx16G "
                "-jar path/to/gatk_4 "
                "-R reference_genome "
                "-T VariantFiltration "
                "--read_filter BadCigar "
                "--read_filter NotPrimaryAlignment "
                "-o tests/assets/jobs/test_dataset/gatk_var_calling/test_user_sample_id_filter_snp.vcf "
                "-l INFO "
                "-U LENIENT_VCF_PROCESSING "
                "-V tests/assets/jobs/test_dataset/gatk_var_calling/test_user_sample_id_raw_snp.vcf.gz "
                "--filterExpression 'QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0' "
                "--filterName 'SNP_FILTER'",
                job_name='var_filtration',
                mem=16,
                working_dir='tests/assets/jobs/test_dataset/gatk_var_calling'
            )
            _test_bgzip_and_tabix(e, 'tests/assets/jobs/test_dataset/gatk_var_calling/test_user_sample_id_filter_snp.vcf')


def _test_bgzip_and_tabix(executor, vcf_file):
    assert executor.call_args_list[1] == call(
        'path/to/bgzip -f ' + vcf_file, cpus=1, job_name='bgzip', mem=8,
        working_dir='tests/assets/jobs/test_dataset/gatk_var_calling'
    )
    assert executor.call_args_list[2] == call(
        'path/to/tabix -f -p vcf ' + vcf_file + '.gz', cpus=1, job_name='tabix', mem=8,
        working_dir='tests/assets/jobs/test_dataset/gatk_var_calling'
    )
