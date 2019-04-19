import os
from os.path import join, dirname, abspath
from unittest.mock import patch

from egcg_core.constants import ELEMENT_LANE, ELEMENT_NB_READS_CLEANED, ELEMENT_RUN_NAME, ELEMENT_RUN_ELEMENT_ID, \
    ELEMENT_PROJECT_ID

from analysis_driver.pipelines.qc_gatk4 import SplitBWA, SplitHaplotypeCaller, SplitFastqStage, FastqIndex, \
    MergeBamAndDup, PostAlignmentScatter
from tests.test_pipelines.test_variant_calling import TestVariantCalling

patch_executor = patch('analysis_driver.pipelines.variant_calling_gatk4.executor.execute')


class TestGATK4(TestVariantCalling):
    def setUp(self):
        super().setUp()
        self.current_wd = os.curdir
        os.chdir(dirname(dirname(__file__)))

    def tearDown(self):
        super().tearDown()
        os.chdir(self.current_wd)


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
            stage._run()
            commands = [
                'gunzip -c fastq_file1.gz | path/to/pbgzip -n 16  -c /dev/stdin > ' + abspath(join(self.assets_path, 'indexed_fastq_file1.gz')),
                'gunzip -c fastq_file2.gz | path/to/pbgzip -n 16  -c /dev/stdin > ' + abspath(join(self.assets_path, 'indexed_fastq_file2.gz'))
            ]
            e.assert_called_with(
                *commands,
                job_name='compress_fastq',
                log_commands=False,
                mem=8,
                working_dir='tests/assets/jobs/test_dataset')


class TestSplitBWA(TestGATK4):

    def test_bwa_command(self):
        stage = SplitBWA(dataset=self.dataset)
        cmd = stage.bwa_command(
            ['file_R1.fastq.gz', 'file_R2.fastq.gz'], 'reference.fa',
            'expected_output_bam', {'ID': '1', 'SM': 'sample1', 'PL': 'illumina'}, (1, 10000))
        exp = 'set -o pipefail; '\
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
                             cpus=2, job_name='bwa_split', mem=12, working_dir='tests/assets/jobs/test_dataset')

        mcommand.assert_any_call(
            chunk=(1, 100000000),
            expected_output_bam=join('tests/assets/jobs/test_dataset/split_alignment/a_run_1_name_sort_1_100000000.bam'),
            fastq_pair=[
                join(self.assets_path, 'indexed_fastq_file1.gz'), join(self.assets_path, 'indexed_fastq_file2.gz')
            ],
            read_group={'ID': 'a_run_1', 'PU': 'a_run_1', 'SM': 'test_user_sample_id', 'PL': 'illumina'},
            reference='reference_genome'
            )
        assert mcommand.call_count == 5


class TestMergeBamAndDup(TestGATK4):

    def _run(self):
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
                      'path/to/bamcat level=0 tmpfile=tests/assets/jobs/test_dataset/cat_tmp/test_dataset ' \
                      '`cat tests/assets/jobs/test_dataset/test_dataset_all_bam_files.list` | ' \
                      'path/to/sortmapdup threads=16 tmpfile=tests/assets/jobs/test_dataset/dup_tmp/test_dataset ' \
                      'indexfilename=tests/assets/jobs/test_dataset/test_dataset.bam.bai ' \
                      '> tests/assets/jobs/test_dataset/test_dataset.bam'
            e.assert_called_with(
                command,
                cpus=6,
                job_name='merge_dup_bam',
                mem=36,
                working_dir='tests/assets/jobs/test_dataset'
            )


class TestPostAlignmentScatter(TestGATK4):

    def test_split_genome_files(self):
        self.dataset.reference_genome = join(self.assets_path, 'genome.fa')
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
            ofile.readlines() == [
                'bigchr2\t40000000\t45000000\n',
                'smchr1\t0\t10000\n',
                'smchr2\t0\t10000\n',
                'smchr3\t0\t10000\n',
                'smchr4\t0\t10000\n'
            ]



class TestScatterBaseRecalibrator(TestGATK4):

    def test_base_recalibrator_cmd(self):
        pass

    def test_run(self):
        pass


class TestGatherBQSRReport(TestGATK4):
    def test_run(self):
        pass


class TestScatterApplyBQSR(TestGATK4):

    def test_apply_bqsr_cmd(self):
        pass

    def test_run(self):
        pass


class TestGatherRecalBam(TestGATK4):

    def test_run(self):
        pass


class TestSplitHaplotypeCaller(TestGATK4):

    def test_haplotype_caller_cmd(self, chunk, region_file):
        pass

    def test_run(self):
        pass


class TestGatherGVCF(TestGATK4):
    def test_run(self):
        pass


class TestSplitGenotypeGVCFs(TestGATK4):

    def test_genotypegvcf_cmd(self):
        pass

    def test_run(self):
        pass


class TestGatherVCF(TestGATK4):

    def test_run(self):
        pass


class TestVariantAnnotation(TestGATK4):

    def test_run(self):
        pass


class TestVQSRFiltrationSNPs(TestGATK4):
    def test_run(self):
        pass


class TestVQSRFiltrationIndels(TestGATK4):
    def test_run(self):
        pass


class TestApplyVQSRSNPs(TestGATK4):
    def test_run(self):
        pass


class TestApplyVQSRIndels(TestGATK4):
    def test_run(self):
        pass


class TestMergeVariants(TestGATK4):

    def test_run(self):
        pass


class TestSelectVariants(TestGATK4):

    def test_run(self):
        pass

class TestSNPsFiltration(TestGATK4):
    def test_run(self):
        pass


class TestIndelsFiltration(TestGATK4):
    def test_run(self):
        pass


class TestSplitHaplotypeCaller(TestGATK4):

    def test_split_genome_in_chunks(self):
        stage = SplitHaplotypeCaller(dataset=self.dataset, bam_file='a_bam_file.bam')
        self.dataset.reference_genome = join(self.assets_path, 'genome.fa')
        print('\n'.join(str(e) for e in stage.split_genome_in_chunks()))

    def test_run(self):
        pass



