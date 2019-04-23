import os
import shutil
from os.path import join, dirname, abspath, isfile, isdir
from unittest.mock import patch, Mock

import pytest
from egcg_core.constants import ELEMENT_LANE, ELEMENT_NB_READS_CLEANED, ELEMENT_RUN_NAME, ELEMENT_RUN_ELEMENT_ID, \
    ELEMENT_PROJECT_ID

from analysis_driver.exceptions import AnalysisDriverError
from analysis_driver.pipelines.human_variant_calling_gatk4 import VQSRFiltrationSNPs, VQSRFiltrationIndels, \
    ApplyVQSRSNPs, ApplyVQSRIndels
from analysis_driver.pipelines.qc_gatk4 import SplitBWA, SplitHaplotypeCaller, SplitFastqStage, FastqIndex, \
    MergeBamAndDup, PostAlignmentScatter, GatherVCF, MergeVariants, SelectSNPs, SelectIndels, IndelsFiltration, \
    SNPsFiltration, GATK4FilePath
from analysis_driver.pipelines.variant_calling_gatk4 import ScatterBaseRecalibrator, PostAlignmentScatterVC, \
    GatherBQSRReport, ScatterApplyBQSR, GatherRecalBam, SplitHaplotypeCallerVC, GatherGVCF, SplitGenotypeGVCFs, \
    VariantAnnotation
from tests.test_pipelines.test_variant_calling import TestVariantCalling


patch_executor = patch('analysis_driver.pipelines.variant_calling_gatk4.executor.execute')


class TestGATK4(TestVariantCalling):
    def setUp(self):
        super().setUp()
        self.dataset.reference_genome = join(self.assets_path, 'genome.fa')
        self.current_wd = os.getcwd()
        os.chdir(dirname(dirname(dirname(__file__))))

    def tearDown(self):
        super().tearDown()
        # cleanup test_dataset job dir
        test_dataset_dir = join(self.assets_path, 'jobs', 'test_dataset')
        for d in os.listdir(test_dataset_dir):
            if isdir(join(test_dataset_dir, d)):
                shutil.rmtree(join(test_dataset_dir, d))
            else:
                os.remove(join(test_dataset_dir, d))
        os.chdir(self.current_wd)


class TestGATK4FilePath(TestGATK4):
    def setUp(self):
        super().setUp()
        self.stage = GATK4FilePath(dataset=self.dataset)

    def test_file_paths(self):
        assert self.stage.hard_filtered_vcf == 'tests/assets/jobs/test_dataset/gatk4/test_user_sample_id_hard_filter.vcf.gz'
        assert self.stage.vqsr_filtered_vcf == 'tests/assets/jobs/test_dataset/gatk4/test_user_sample_id_vqsr.vcf.gz'
        assert self.stage.vqsr_datasets == {
            'hapmap': '/path/to/hapmap_annotation',
            'omni': '/path/to/omni_annotation',
            'thousand_genomes': '/path/to/1000g_annotation',
            'dbsnp': '/path/to/dbsnp_annotation',
            'mills': '/path/to/mills_annotation'
        }
        with pytest.raises(AnalysisDriverError):
            self.stage.dataset.genome_version = 'do_not_exist'
            print(self.stage.vqsr_datasets)


class TestSplitFastqStage(TestGATK4):

    def test_chunks_from_fastq(self):
        stage = SplitFastqStage(dataset=self.dataset)
        chunks = stage.chunks_from_fastq([os.path.join(self.assets_path, 'indexed_fastq_file1.gz')])
        assert chunks == [
            (1, 100000000), (100000001, 200000000), (200000001, 300000000), (300000001, 400000000),
            (400000001, 500000000)
        ]


class TestFastqIndex(TestGATK4):

    def test_run(self):
        self.dataset.run_elements = [
            {ELEMENT_RUN_ELEMENT_ID: 'a_run_1', ELEMENT_PROJECT_ID: 'a_project', ELEMENT_NB_READS_CLEANED: 1,
             ELEMENT_RUN_NAME: 'a_run', ELEMENT_LANE: 1}
        ]
        index_fastq_files = [join(self.assets_path, 'indexed_fastq_file1.gz'),
                             join(self.assets_path, 'indexed_fastq_file2.gz')]
        fastq_files = ['fastq_file1.gz', 'fastq_file2.gz']
        stage = FastqIndex(dataset=self.dataset)
        with patch.object(FastqIndex, '_find_fastqs_for_run_element', return_value=fastq_files), \
                patch.object(FastqIndex, '_indexed_fastq_for_run_element', return_value=index_fastq_files), \
                patch_executor as e:
            e.return_value = Mock(join=Mock(return_value=0))
            stage._run()
            commands = [
                'gunzip -c fastq_file1.gz | path/to/pbgzip -n 16  -c /dev/stdin > ' + abspath(
                    join(self.assets_path, 'indexed_fastq_file1.gz')),
                'gunzip -c fastq_file2.gz | path/to/pbgzip -n 16  -c /dev/stdin > ' + abspath(
                    join(self.assets_path, 'indexed_fastq_file2.gz'))
            ]
            e.assert_any_call(
                *commands,
                job_name='compress_fastq',
                log_commands=False,
                mem=8,
                working_dir='tests/assets/jobs/test_dataset/slurm_and_logs'
            )
            commands = [
                'path/to/grabix index ' + join(self.assets_path, 'indexed_fastq_file1.gz'),
                'path/to/grabix index ' + join(self.assets_path, 'indexed_fastq_file2.gz')
            ]
            e.assert_any_call(
                *commands, job_name='index_fastq', working_dir='tests/assets/jobs/test_dataset/slurm_and_logs'
            )



class TestSplitBWA(TestGATK4):

    def test_bwa_command(self):
        stage = SplitBWA(dataset=self.dataset)
        cmd = stage.bwa_command(
            ['file_R1.fastq.gz', 'file_R2.fastq.gz'], 'reference.fa',
            'expected_output_bam', {'ID': '1', 'SM': 'sample1', 'PL': 'illumina'}, (1, 10000))
        exp = 'set -o pipefail; ' \
              r'''path/to/bwa_1.1 mem -K 100000000 -Y -R '@RG\tID:1\tPL:illumina\tSM:sample1' -M -t 2 reference.fa ''' \
              '<(path/to/grabix grab file_R1.fastq.gz 1 10000) <(path/to/grabix grab file_R2.fastq.gz 1 10000) | ' \
              'path/to/samtools_1.3.1 sort -n -m 1G -O bam -T expected_output_bam -o expected_output_bam -'
        assert cmd == exp

    def test_run(self):
        self.dataset.run_elements = [
            {ELEMENT_RUN_ELEMENT_ID: 'a_run_1', ELEMENT_NB_READS_CLEANED: 1, ELEMENT_RUN_NAME: 'a_run', ELEMENT_LANE: 1}
        ]
        index_fastq_files = [join(self.assets_path, 'indexed_fastq_file1.gz'),
                             join(self.assets_path, 'indexed_fastq_file2.gz')]

        stage = SplitBWA(dataset=self.dataset)
        with patch.object(SplitBWA, '_indexed_fastq_for_run_element', return_value=index_fastq_files), \
                patch_executor as e, patch.object(SplitBWA, 'bwa_command', return_value='command_bwa') as mcommand:
            stage._run()
        e.assert_called_with('command_bwa', 'command_bwa', 'command_bwa', 'command_bwa', 'command_bwa',
                             cpus=2, job_name='bwa_split', mem=12,
                             working_dir='tests/assets/jobs/test_dataset/slurm_and_logs')

        mcommand.assert_any_call(
            chunk=(1, 100000000),
            expected_output_bam=join(
                'tests/assets/jobs/test_dataset/split_alignment/a_run_1_name_sort_1_100000000.bam'),
            fastq_pair=[
                join(self.assets_path, 'indexed_fastq_file1.gz'), join(self.assets_path, 'indexed_fastq_file2.gz')
            ],
            read_group={'ID': 'a_run_1', 'PU': 'a_run_1', 'SM': 'test_user_sample_id', 'PL': 'illumina'},
            reference=self.dataset.reference_genome
        )
        assert mcommand.call_count == 5


class TestMergeBamAndDup(TestGATK4):

    def test_run(self):
        self.dataset.run_elements = [
            {ELEMENT_RUN_ELEMENT_ID: 'a_run_1', ELEMENT_NB_READS_CLEANED: 1, ELEMENT_RUN_NAME: 'a_run', ELEMENT_LANE: 1}
        ]
        index_fastq_files = [join(self.assets_path, 'indexed_fastq_file1.gz'),
                             join(self.assets_path, 'indexed_fastq_file2.gz')]
        stage = MergeBamAndDup(dataset=self.dataset)

        with patch.object(MergeBamAndDup, '_indexed_fastq_for_run_element', return_value=index_fastq_files), \
                patch_executor as e:
            stage._run()
            command = 'set -o pipefail; ' \
                      'path/to/bamcat level=0 tmpfile=%s ' \
                      '`cat tests/assets/jobs/test_dataset/test_dataset_all_bam_files.list` | ' \
                      'path/to/sortmapdup threads=16 tmpfile=%s ' \
                      'indexfilename=tests/assets/jobs/test_dataset/test_dataset.bam.bai ' \
                      '> tests/assets/jobs/test_dataset/test_dataset.bam' % (
                          os.path.join(stage.dir_to_delete[0], 'test_dataset'),
                          os.path.join(stage.dir_to_delete[1], 'test_dataset')
                      )
            e.assert_called_with(
                command,
                cpus=6,
                job_name='merge_dup_bam',
                mem=36,
                working_dir='tests/assets/jobs/test_dataset/slurm_and_logs'
            )


class TestPostAlignmentScatter(TestGATK4):

    def test_split_genome_in_chunks1(self):
        stage = PostAlignmentScatter(dataset=self.dataset)
        exp_chunks = [
            [('bigchr1', 0, 15000000)],
            [('bigchr2', 0, 20000000)],
            [('bigchr2', 20000000, 40000000)],
            [
                ('bigchr2', 40000000, 45000000), ('smchr1', 0, 10000), ('smchr2', 0, 10000),
                ('smchr3', 0, 10000), ('smchr4', 0, 10000)
            ],
            [('bigchr3', 0, 20000000)],
            [('bigchr3', 20000000, 40000000)],
            [('bigchr3', 40000000, 60000000)],
            [('bigchr3', 60000000, 76000000)]
        ]
        assert exp_chunks == stage.split_genome_in_chunks()

    def test_split_genome_in_chunks2(self):
        ref_file = join(self.assets_path, 'genome_many_contigs.fa')
        with open(ref_file + '.fai', 'w') as open_file:
            for i in range(2400):
                open_file.write('contig_%s\t1000\t-\t-\t-\n' % i)

        self.dataset.reference_genome = ref_file
        stage = PostAlignmentScatter(dataset=self.dataset)

        obs_chunks = stage.split_genome_in_chunks()
        # 3 sets of chunks because limit of 1000 chunk per file
        assert len(obs_chunks) == 3
        os.remove(ref_file + '.fai')

    def test_split_genome_files(self):
        stage = PostAlignmentScatter(dataset=self.dataset)
        split_dir = 'tests/assets/jobs/test_dataset/post_alignment_split'
        expected_output = {
            ('bigchr1', 0, 15000000): join(split_dir, 'test_dataset_region_bigchr1-0-15000000.bed'),
            ('bigchr2', 0, 20000000): join(split_dir, 'test_dataset_region_bigchr2-0-20000000.bed'),
            ('bigchr2', 20000000, 40000000): join(split_dir, 'test_dataset_region_bigchr2-20000000-40000000.bed'),
            ('bigchr2', 40000000, 45000000): join(split_dir, 'test_dataset_region_bigchr2-40000000-45000000.bed'),
            ('bigchr3', 0, 20000000): join(split_dir, 'test_dataset_region_bigchr3-0-20000000.bed'),
            ('bigchr3', 20000000, 40000000): join(split_dir, 'test_dataset_region_bigchr3-20000000-40000000.bed'),
            ('bigchr3', 40000000, 60000000): join(split_dir, 'test_dataset_region_bigchr3-40000000-60000000.bed'),
            ('bigchr3', 60000000, 76000000): join(split_dir, 'test_dataset_region_bigchr3-60000000-76000000.bed')
        }

        assert stage.split_genome_files() == expected_output

        # Small chroms are added along the larger ones in order
        with open(join(split_dir, 'test_dataset_region_bigchr2-40000000-45000000.bed')) as ofile:
            assert ofile.readlines() == [
                'bigchr2\t40000000\t45000000\n',
                'smchr1\t0\t10000\n',
                'smchr2\t0\t10000\n',
                'smchr3\t0\t10000\n',
                'smchr4\t0\t10000\n'
            ]


class TestPostAlignmentScatterVC(TestGATK4):

    def test_split_genome_chromosomes(self):
        stage = PostAlignmentScatterVC(dataset=self.dataset)
        assert stage.split_genome_chromosomes() == [
            ['bigchr1', 'bigchr2', 'smchr1', 'smchr2', 'smchr3', 'smchr4'],
            ['bigchr3']
        ]


class TestScatterBaseRecalibrator(TestGATK4):

    def test_base_recalibrator_cmd(self):
        stage = ScatterBaseRecalibrator(dataset=self.dataset)
        obs_cmd = stage.base_recalibrator_cmd(['chr1'])
        cmd = 'path/to/gatk ' \
              '--java-options "-Djava.io.tmpdir=%s -XX:+UseSerialGC -Xmx6G" BaseRecalibrator ' \
              '--output tests/assets/jobs/test_dataset/post_alignment_split/test_dataset_base_recal_grp_chr1.grp ' \
              '--input tests/assets/jobs/test_dataset/test_dataset.bam  ' \
              '--reference %s --known-sites  /path/to/dbsnp.vcf.gz ' \
              '--intervals chr1' % (stage.dir_to_delete[0], self.dataset.reference_genome)
        assert obs_cmd == cmd

    def test_run(self):
        stage = ScatterBaseRecalibrator(dataset=self.dataset)
        with patch.object(ScatterBaseRecalibrator, 'base_recalibrator_cmd', return_value='recal_cmd') as mcommand, \
                patch_executor as e:
            stage._run()
            e.assert_called_with(
                'recal_cmd', 'recal_cmd',
                cpus=1, job_name='gatk_base_recal', mem=6, working_dir='tests/assets/jobs/test_dataset/slurm_and_logs'
            )
            assert mcommand.call_count == 2


class TestGatherBQSRReport(TestGATK4):

    def test_run(self):
        stage = GatherBQSRReport(dataset=self.dataset)
        with patch_executor as e:
            stage._run()
            e.assert_called_with(
                'path/to/gatk --java-options "-Djava.io.tmpdir=%s -XX:+UseSerialGC -Xmx6G" GatherBQSRReports '
                '--output tests/assets/jobs/test_dataset/gatk4/test_user_sample_id.grp '
                '--input tests/assets/jobs/test_dataset/post_alignment_split/test_dataset_bqsr_reports.list'
                ' ' % (stage.dir_to_delete[0]), cpus=1, job_name='gather_bqsr', mem=6,
                working_dir='tests/assets/jobs/test_dataset/slurm_and_logs')


class TestScatterApplyBQSR(TestGATK4):

    def test_apply_bqsr_cmd(self):
        stage = ScatterApplyBQSR(dataset=self.dataset)
        obs_cmd = stage.apply_bqsr_cmd(['chr1'])
        exp_cmd = 'path/to/gatk --java-options "-Djava.io.tmpdir=%s -XX:+UseSerialGC -Xmx6G" ApplyBQSR ' \
                  '--output tests/assets/jobs/test_dataset/post_alignment_split/test_dataset_recal_chr1.bam ' \
                  '--input tests/assets/jobs/test_dataset/test_dataset.bam  ' \
                  '--reference %s ' \
                  '--bqsr-recal-file tests/assets/jobs/test_dataset/gatk4/test_user_sample_id.grp ' \
                  '--jdk-deflater --jdk-inflater ' \
                  '--static-quantized-quals 10 --static-quantized-quals 20 --static-quantized-quals 30 ' \
                  '--static-quantized-quals 40 --intervals chr1' % (
                      stage.dir_to_delete[0], self.dataset.reference_genome
                  )
        assert obs_cmd == exp_cmd

    def test_run(self):
        stage = ScatterApplyBQSR(dataset=self.dataset)
        with patch.object(ScatterApplyBQSR, 'apply_bqsr_cmd', return_value='bqsr_cmd') as mcommand, \
                patch_executor as e:
            stage._run()
            e.assert_called_with(
                'bqsr_cmd', 'bqsr_cmd', 'bqsr_cmd',
                cpus=1, job_name='apply_bqsr', mem=6, working_dir='tests/assets/jobs/test_dataset/slurm_and_logs'
            )
            # The two sets of chrom and the unmapped
            mcommand.assert_any_call(['bigchr1', 'bigchr2', 'smchr1', 'smchr2', 'smchr3', 'smchr4'])
            mcommand.assert_any_call(['bigchr3'])
            mcommand.assert_any_call(['unmapped'])
            assert mcommand.call_count == 3


class TestGatherRecalBam(TestGATK4):

    def test_run(self):
        stage = GatherRecalBam(dataset=self.dataset)
        with patch_executor as e:
            stage._run()
        exp_cmd = 'path/to/java_8 -Djava.io.tmpdir=%s -XX:+UseSerialGC -Xmx6G -jar path/to/picard GatherBamFiles ' \
                  'INPUT=tests/assets/jobs/test_dataset/post_alignment_split/test_dataset_recal_bam.list ' \
                  'OUTPUT=tests/assets/jobs/test_dataset/gatk4/test_user_sample_id_recal.bam ' \
                  'VALIDATION_STRINGENCY=LENIENT CREATE_INDEX=true' % stage.dir_to_delete[0]
        e.assert_called_with(exp_cmd, cpus=1, job_name='gather_recal_bam', mem=6,
                             working_dir='tests/assets/jobs/test_dataset/slurm_and_logs')


class TestSplitHaplotypeCaller(TestGATK4):

    def test_haplotype_caller_cmd(self):
        stage = SplitHaplotypeCaller(dataset=self.dataset, bam_file='a_bam_file.bam')
        obs_cmd = stage.haplotype_caller_cmd(('chr1', 1, 1000), 'path/to/region/file')
        exp_cmd = 'path/to/gatk --java-options "-Djava.io.tmpdir=%s ' \
                  '-XX:+UseSerialGC -Xmx6G" HaplotypeCaller ' \
                  '--output tests/assets/jobs/test_dataset/post_alignment_split/test_dataset_haplotype_caller_chr1-1-1000.vcf.gz ' \
                  '--input a_bam_file.bam  ' \
                  '--reference %s  ' \
                  '--sample-ploidy 2 --intervals path/to/region/file --annotation BaseQualityRankSumTest ' \
                  '--annotation ClippingRankSumTest --annotation Coverage --annotation DepthPerAlleleBySample ' \
                  '--annotation DepthPerSampleHC --annotation FisherStrand --annotation MappingQuality ' \
                  '--annotation MappingQualityRankSumTest --annotation MappingQualityZero --annotation QualByDepth ' \
                  '--annotation ReadPosRankSumTest --annotation RMSMappingQuality --dbsnp /path/to/dbsnp.vcf.gz'
        assert exp_cmd % (stage.dir_to_delete[0], self.dataset.reference_genome) == obs_cmd

    def test_run(self):
        stage = SplitHaplotypeCaller(dataset=self.dataset, bam_file='a_bam_file.bam')
        with patch.object(SplitHaplotypeCaller, 'haplotype_caller_cmd', return_value='hp_cmd') as mcommand, \
                patch_executor as e:
            stage._run()
            e.assert_called_with(*['hp_cmd'] * 8, cpus=1, job_name='split_hap_call', mem=6,
                                 working_dir='tests/assets/jobs/test_dataset/slurm_and_logs')
            mcommand.assert_any_call(
                ('bigchr1', 0, 15000000),
                'tests/assets/jobs/test_dataset/post_alignment_split/test_dataset_region_bigchr1-0-15000000.bed'
            )
            assert mcommand.call_count == 8


class TestSplitHaplotypeCallerVC(TestGATK4):
    def test_haplotype_caller_cmd(self):
        stage = SplitHaplotypeCallerVC(dataset=self.dataset, bam_file='a_bam_file.bam')
        obs_cmd = stage.haplotype_caller_cmd(('chr1', 1, 1000), 'path/to/region/file')
        exp_cmd = 'path/to/gatk --java-options "-Djava.io.tmpdir=%s ' \
                  '-XX:+UseSerialGC -Xmx6G" HaplotypeCaller ' \
                  '--output tests/assets/jobs/test_dataset/post_alignment_split/test_dataset_haplotype_caller_chr1-1-1000.g.vcf.gz ' \
                  '--input a_bam_file.bam  ' \
                  '--reference %s  ' \
                  '--sample-ploidy 2 --emit-ref-confidence GVCF --intervals path/to/region/file ' \
                  '--annotation BaseQualityRankSumTest --annotation ClippingRankSumTest --annotation Coverage ' \
                  '--annotation DepthPerAlleleBySample --annotation DepthPerSampleHC --annotation FisherStrand ' \
                  '--annotation MappingQuality --annotation MappingQualityRankSumTest --annotation MappingQualityZero '\
                  '--annotation QualByDepth --annotation ReadPosRankSumTest --annotation RMSMappingQuality ' \
                  '--dbsnp /path/to/dbsnp.vcf.gz'
        assert exp_cmd % (stage.dir_to_delete[0], self.dataset.reference_genome) == obs_cmd


class TestGatherGVCF(TestGATK4):
    def test_run(self):
        stage = GatherGVCF(dataset=self.dataset)
        with patch_executor as e:
            e.return_value = Mock(join=Mock(return_value=0))
            stage._run()
            exp_command = 'path/to/gatk --java-options "-Djava.io.tmpdir=%s -XX:+UseSerialGC -Xmx6G" GatherVcfs ' \
                          '--OUTPUT tests/assets/jobs/test_dataset/gatk4/test_user_sample_id.g.vcf.gz ' \
                          '--INPUT tests/assets/jobs/test_dataset/post_alignment_split/test_dataset_g.vcf.list '
            e.assert_any_call(exp_command % stage.dir_to_delete[0], cpus=1, job_name='gather_hap_call', mem=6,
                              working_dir='tests/assets/jobs/test_dataset/slurm_and_logs')
            e.assert_any_call('path/to/tabix -f -p vcf tests/assets/jobs/test_dataset/gatk4/test_user_sample_id.g.vcf.gz',
                              cpus=1, job_name='tabix', mem=8, working_dir='tests/assets/jobs/test_dataset/slurm_and_logs')
            assert e.call_count == 2


class TestSplitGenotypeGVCFs(TestGATK4):

    def test_genotypegvcf_cmd(self):
        stage = SplitGenotypeGVCFs(dataset=self.dataset)
        obs_cmd = stage.genotypegvcf_cmd(('chr1', 1, 1000), 'path/to/region/file')
        exp_cmd = 'path/to/gatk --java-options "-Djava.io.tmpdir=%s -XX:+UseSerialGC -Xmx6G" GenotypeGVCFs ' \
                  '--output tests/assets/jobs/test_dataset/post_alignment_split/test_dataset_genotype_gvcf_chr1-1-1000.vcf.gz  ' \
                  '--reference %s ' \
                  '--variant tests/assets/jobs/test_dataset/post_alignment_split/test_dataset_haplotype_caller_chr1-1-1000.g.vcf.gz ' \
                  '--sample-ploidy 2 --intervals path/to/region/file --annotation BaseQualityRankSumTest ' \
                  '--annotation ClippingRankSumTest --annotation Coverage --annotation DepthPerAlleleBySample ' \
                  '--annotation DepthPerSampleHC --annotation FisherStrand --annotation MappingQuality ' \
                  '--annotation MappingQualityRankSumTest --annotation MappingQualityZero --annotation QualByDepth ' \
                  '--annotation ReadPosRankSumTest --annotation RMSMappingQuality'
        assert exp_cmd % (stage.dir_to_delete[0], self.dataset.reference_genome) == obs_cmd

    def test_run(self):
        stage = SplitGenotypeGVCFs(dataset=self.dataset)
        with patch.object(SplitGenotypeGVCFs, 'genotypegvcf_cmd', return_value='gg_cmd') as mcommand, \
                patch_executor as e:
            stage._run()
            e.assert_called_with(*['gg_cmd'] * 8, cpus=1, job_name='split_genotype_call', mem=6,
                                 working_dir='tests/assets/jobs/test_dataset/slurm_and_logs')
            mcommand.assert_any_call(
                ('bigchr1', 0, 15000000),
                'tests/assets/jobs/test_dataset/post_alignment_split/test_dataset_region_bigchr1-0-15000000.bed'
            )
            assert mcommand.call_count == 8


class TestGatherVCF(TestGATK4):

    def test_run(self):
        stage = GatherVCF(dataset=self.dataset)
        with patch_executor as e:
            e.return_value = Mock(join=Mock(return_value=0))
            stage._run()
            exp_cmd = 'path/to/gatk --java-options "-Djava.io.tmpdir=%s -XX:+UseSerialGC -Xmx6G" GatherVcfs ' \
                      '--OUTPUT tests/assets/jobs/test_dataset/gatk4/test_user_sample_id.vcf.gz ' \
                      '--INPUT tests/assets/jobs/test_dataset/post_alignment_split/test_dataset_vcf.list ' % stage.dir_to_delete[0]
            e.assert_any_call(exp_cmd, cpus=1, job_name='gather_geno_call', mem=6,
                              working_dir='tests/assets/jobs/test_dataset/slurm_and_logs')


class TestVariantAnnotation(TestGATK4):

    def test_run(self):
        stage = VariantAnnotation(dataset=self.dataset)
        with patch_executor as e:
            e.return_value = Mock(join=Mock(return_value=0))
            stage._run()
            exp_cmd = 'set -o pipefail; path/to/java_8 -Djava.io.tmpdir=%s -XX:+UseSerialGC -Xmx20G ' \
                      '-jar path/to/snpEff eff -hgvs -noLog -i vcf -o vcf ' \
                      '-csvStats tests/assets/jobs/test_dataset/gatk4/test_user_sample_id_snpseff.csv ' \
                      '-s tests/assets/jobs/test_dataset/gatk4/test_user_sample_id_snpseff.html ' \
                      'snpEffdb tests/assets/jobs/test_dataset/gatk4/test_user_sample_id.vcf.gz | ' \
                      'path/to/bgzip --threads 16 -c > ' \
                      'tests/assets/jobs/test_dataset/gatk4/test_user_sample_id_snpseff.vcf.gz'
            e.assert_any_call(exp_cmd % stage.dir_to_delete[0], cpus=2, job_name='snpeff', mem=20,
                              working_dir='tests/assets/jobs/test_dataset/slurm_and_logs')
            e.assert_any_call('path/to/tabix -f -p vcf '
                              'tests/assets/jobs/test_dataset/gatk4/test_user_sample_id_snpseff.vcf.gz',
                              cpus=1, job_name='tabix', mem=8,
                              working_dir='tests/assets/jobs/test_dataset/slurm_and_logs')
            assert e.call_count == 2


class TestVQSRFiltrationSNPs(TestGATK4):
    def test_run(self):
        stage = VQSRFiltrationSNPs(dataset=self.dataset)
        with patch_executor as e:
            stage._run()
            exp_cmd = 'path/to/gatk --java-options "-Djava.io.tmpdir=%s -XX:+UseSerialGC -Xmx16G" VariantRecalibrator ' \
                      '--output tests/assets/jobs/test_dataset/gatk4/test_user_sample_id_vqsr_snps_recall.vcf.gz  ' \
                      '--reference %s ' \
                      '-V tests/assets/jobs/test_dataset/gatk4/test_user_sample_id_snpseff.vcf.gz ' \
                      '--resource:hapmap,known=false,training=true,truth=true,prior=15.0 /path/to/hapmap_annotation ' \
                      '--resource:omni,known=false,training=true,truth=false,prior=12.0 /path/to/omni_annotation ' \
                      '--resource:1000G,known=false,training=true,truth=false,prior=10.0 /path/to/1000g_annotation ' \
                      '--resource:dbsnp,known=true,training=false,truth=false,prior=2.0 /path/to/dbsnp_annotation ' \
                      '--max-gaussians 4 -tranche 100.0 -tranche 99.99 -tranche 99.98 -tranche 99.97 -tranche 99.96 ' \
                      '-tranche 99.95 -tranche 99.94 -tranche 99.93 -tranche 99.92 -tranche 99.91 -tranche 99.9 ' \
                      '-tranche 99.8 -tranche 99.7 -tranche 99.6 -tranche 99.5 -tranche 99.0 -tranche 98.0 ' \
                      '-tranche 90.0 -an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR -mode SNP ' \
                      '--tranches-file tests/assets/jobs/test_dataset/gatk4/test_user_sample_id_vqsr_snps.tranches ' \
                      '--rscript-file tests/assets/jobs/test_dataset/gatk4/test_user_sample_id_vqsr_snps.R'
            e.assert_called_with(exp_cmd % (stage.dir_to_delete[0], self.dataset.reference_genome), cpus=1,
                                 job_name='vqsr_snp', mem=16,
                                 working_dir='tests/assets/jobs/test_dataset/slurm_and_logs')


class TestVQSRFiltrationIndels(TestGATK4):
    def test_run(self):
        stage = VQSRFiltrationIndels(dataset=self.dataset)
        with patch_executor as e:
            stage._run()
            exp_cmd = 'path/to/gatk --java-options "-Djava.io.tmpdir=%s -XX:+UseSerialGC -Xmx16G" VariantRecalibrator ' \
                      '--output tests/assets/jobs/test_dataset/gatk4/test_user_sample_id_vqsr_indels_recall.vcf.gz  ' \
                      '--reference %s -V tests/assets/jobs/test_dataset/gatk4/test_user_sample_id_snpseff.vcf.gz ' \
                      '--resource:mills,known=false,training=true,truth=true,prior=12.0 /path/to/mills_annotation ' \
                      '--resource:dbsnp,known=true,training=false,truth=false,prior=2.0 /path/to/dbsnp_annotation ' \
                      '--max-gaussians 4 -tranche 100.0 -tranche 99.99 -tranche 99.98 -tranche 99.97 -tranche 99.96 ' \
                      '-tranche 99.95 -tranche 99.94 -tranche 99.93 -tranche 99.92 -tranche 99.91 -tranche 99.9 ' \
                      '-tranche 99.8 -tranche 99.7 -tranche 99.6 -tranche 99.5 -tranche 99.0 -tranche 98.0 ' \
                      '-tranche 90.0 -an QD -an MQRankSum -an ReadPosRankSum -an FS -an SOR -an DP -mode INDEL ' \
                      '--tranches-file tests/assets/jobs/test_dataset/gatk4/test_user_sample_id_vqsr_indels_tranches ' \
                      '--rscript-file tests/assets/jobs/test_dataset/gatk4/test_user_sample_id_vqsr_indels.R'
            e.assert_called_with(exp_cmd % (stage.dir_to_delete[0], self.dataset.reference_genome), cpus=1,
                                 job_name='vqsr_indel', mem=16,
                                 working_dir='tests/assets/jobs/test_dataset/slurm_and_logs')


class TestApplyVQSRSNPs(TestGATK4):
    def test_run(self):
        stage = ApplyVQSRSNPs(dataset=self.dataset)
        with patch_executor as e:
            stage._run()
            exp_cmd = 'path/to/gatk --java-options "-Djava.io.tmpdir=%s -XX:+UseSerialGC -Xmx16G" ApplyVQSR ' \
                      '--output tests/assets/jobs/test_dataset/gatk4/test_user_sample_id_vqsr_snps.vcf.gz  ' \
                      '--reference %s -V tests/assets/jobs/test_dataset/gatk4/test_user_sample_id_snpseff.vcf.gz ' \
                      '-mode SNP ' \
                      '--tranches-file tests/assets/jobs/test_dataset/gatk4/test_user_sample_id_vqsr_snps.tranches ' \
                      '--truth-sensitivity-filter-level 99.0 ' \
                      '--recal-file tests/assets/jobs/test_dataset/gatk4/test_user_sample_id_vqsr_snps_recall.vcf.gz'
            e.assert_called_with(exp_cmd % (stage.dir_to_delete[0], self.dataset.reference_genome), cpus=1,
                                 job_name='apply_vqsr_snps', mem=16,
                                 working_dir='tests/assets/jobs/test_dataset/slurm_and_logs')


class TestApplyVQSRIndels(TestGATK4):
    def test_run(self):
        stage = ApplyVQSRIndels(dataset=self.dataset)
        with patch_executor as e:
            stage._run()
            exp_cmd = 'path/to/gatk --java-options "-Djava.io.tmpdir=%s -XX:+UseSerialGC -Xmx8G" ApplyVQSR ' \
                      '--output tests/assets/jobs/test_dataset/gatk4/test_user_sample_id_vqsr_indels.vcf.gz  ' \
                      '--reference %s ' \
                      '-V tests/assets/jobs/test_dataset/gatk4/test_user_sample_id_snpseff.vcf.gz -mode INDEL ' \
                      '--tranches-file tests/assets/jobs/test_dataset/gatk4/test_user_sample_id_vqsr_indels_tranches ' \
                      '--truth-sensitivity-filter-level 99.0 ' \
                      '--recal-file tests/assets/jobs/test_dataset/gatk4/test_user_sample_id_vqsr_indels_recall.vcf.gz'
            e.assert_called_with(exp_cmd % (stage.dir_to_delete[0], self.dataset.reference_genome), cpus=1,
                                 job_name='apply_vqsr_indels', mem=16,
                                 working_dir='tests/assets/jobs/test_dataset/slurm_and_logs')


class TestSelectSNPs(TestGATK4):
    def test_run(self):
        stage = SelectSNPs(dataset=self.dataset, input_vcf='input_file.vcf.gz')
        with patch_executor as e:
            stage._run()
            exp_cmd = 'path/to/gatk --java-options "-Djava.io.tmpdir=%s -XX:+UseSerialGC -Xmx2G" SelectVariants ' \
                      '--output tests/assets/jobs/test_dataset/gatk4/test_user_sample_id_raw_snp.vcf  ' \
                      '--reference %s ' \
                      '-V input_file.vcf.gz -select-type SNP '
            e.assert_called_with(exp_cmd % (stage.dir_to_delete[0], self.dataset.reference_genome),
                                 job_name='snp_select', mem=16,
                                 working_dir='tests/assets/jobs/test_dataset/slurm_and_logs')


class TestSelectIndels(TestGATK4):
    def test_run(self):
        stage = SelectIndels(dataset=self.dataset, input_vcf='input_file.vcf.gz')
        with patch_executor as e:
            stage._run()
            exp_cmd = 'path/to/gatk --java-options "-Djava.io.tmpdir=%s -XX:+UseSerialGC -Xmx2G" SelectVariants ' \
                      '--output tests/assets/jobs/test_dataset/gatk4/test_user_sample_id_raw_indel.vcf  ' \
                      '--reference %s ' \
                      '-V input_file.vcf.gz -select-type INDEL '
            e.assert_called_with(exp_cmd % (stage.dir_to_delete[0], self.dataset.reference_genome),
                                 job_name='indel_select', mem=16,
                                 working_dir='tests/assets/jobs/test_dataset/slurm_and_logs')


class TestSNPsFiltration(TestGATK4):
    def test_run(self):
        stage = SNPsFiltration(dataset=self.dataset)
        with patch_executor as e:
            stage._run()
            exp_cmd = 'path/to/gatk --java-options "-Djava.io.tmpdir=%s -XX:+UseSerialGC -Xmx2G" VariantFiltration ' \
                      '--output tests/assets/jobs/test_dataset/gatk4/test_user_sample_id_hard_snps.vcf.gz  ' \
                      '--reference %s ' \
                      '-V tests/assets/jobs/test_dataset/gatk4/test_user_sample_id_raw_snp.vcf ' \
                      '--filter-expression \'QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0 || SOR > 3.0\' ' \
                      '--filter-name \'SNP_HARD_FILTER\''
            e.assert_called_with(exp_cmd % (stage.dir_to_delete[0], self.dataset.reference_genome),
                                 job_name='snps_filtration', mem=8,
                                 working_dir='tests/assets/jobs/test_dataset/slurm_and_logs')


class TestIndelsFiltration(TestGATK4):
    def test_run(self):
        stage = IndelsFiltration(dataset=self.dataset)
        with patch_executor as e:
            stage._run()
            exp_cmd = 'path/to/gatk --java-options "-Djava.io.tmpdir=%s -XX:+UseSerialGC -Xmx2G" VariantFiltration ' \
                      '--output tests/assets/jobs/test_dataset/gatk4/test_user_sample_id_hard_filter_indels.vcf.gz  ' \
                      '--reference %s ' \
                      '-V tests/assets/jobs/test_dataset/gatk4/test_user_sample_id_raw_indel.vcf ' \
                      '--filter-expression \'QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0 || SOR > 10.0\' ' \
                      '--filter-name \'INDEL_HARD_FILTER\''
            e.assert_called_with(exp_cmd% (stage.dir_to_delete[0], self.dataset.reference_genome),
                                 job_name='indel_filtration', mem=16,
                                 working_dir='tests/assets/jobs/test_dataset/slurm_and_logs')


class TestMergeVariants(TestGATK4):

    def test_run(self):
        stage = MergeVariants(dataset=self.dataset, vcf_files=['file1.vcf.gz', 'file1.vcf.gz'],
                              output_vcf_file='output.vcf.gz')
        with patch_executor as e:
            stage._run()
            exp_cmd = 'path/to/gatk --java-options "-Djava.io.tmpdir=%s -XX:+UseSerialGC -Xmx2G" MergeVcfs ' \
                      '--OUTPUT output.vcf.gz ' \
                      '--INPUT tests/assets/jobs/test_dataset/file1.vcf.list '
            e.assert_called_with(exp_cmd % stage.dir_to_delete[0], cpus=1, job_name='merge_vqsr', mem=8,
                                 working_dir='tests/assets/jobs/test_dataset/slurm_and_logs')
