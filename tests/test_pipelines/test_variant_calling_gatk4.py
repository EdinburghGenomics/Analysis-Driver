import os
from os.path import join, dirname
from unittest.mock import patch

from egcg_core.constants import ELEMENT_LANE, ELEMENT_NB_READS_CLEANED, ELEMENT_RUN_NAME, ELEMENT_RUN_ELEMENT_ID

from analysis_driver.pipelines.variant_calling_gatk4 import SplitBWA, SplitHaplotypeCaller, SplitFastqStage

from tests.test_pipelines.test_variant_calling import TestVariantCalling

patch_executor = patch('analysis_driver.pipelines.variant_calling_gatk4.executor.execute')


class TestGATK4(TestVariantCalling):
    def setUp(self):
        self.current_wd = os.curdir
        os.chdir(dirname(dirname(__file__)))

    def tearDown(self):
        os.chdir(self.current_wd)


class TestSplitFastqStage(TestGATK4):

    def test_chunks_from_fastq(self):
        stage = SplitFastqStage(dataset=self.dataset)
        chunks = stage.chunks_from_fastq([os.path.join(self.assets_path, 'indexed_fastq_file1.gz')])
        assert chunks == [
            (1, 100000000), (100000001, 200000000), (200000001, 300000000), (300000001, 400000000),
            (400000001, 500000000)
        ]


class TestFastqIndex(TestVariantCalling):

    def test_run(self):
        self.dataset.run_elements = [
            {ELEMENT_RUN_ELEMENT_ID: 'a_run_1', ELEMENT_NB_READS_CLEANED: 1, ELEMENT_RUN_NAME: 'a_run', ELEMENT_LANE: 1}
        ]
        index_fastq_files = [join(self.assets_path, 'indexed_fastq_file1.gz'),
                             join(self.assets_path, 'indexed_fastq_file2.gz')]
        stage = FastqIndex(dataset=self.dataset)
        with patch.object(FastqIndex, '_indexed_fastq_for_run_element', return_value=index_fastq_files), \
             patch_executor as e, patch.object(SplitBWA, 'bwa_command', return_value='command_bwa') as mcommand:
            self.stage._run()


class TestSplitBWA(TestVariantCalling):

    def test_bwa_command(self):
        stage = SplitBWA(dataset=self.dataset)
        cmd = stage.bwa_command(
            ['file_R1.fastq.gz', 'file_R2.fastq.gz'], 'reference.fa',
            'expected_output_bam', {'ID': '1', 'SM': 'sample1', 'PL': 'illumina'}, (1, 10000))
        exp = 'set -o pipefail; '\
              r'''path/to/bwa_1.1 -R '@RG\tID:1\tPL:illumina\tSM:sample1' mem -M  reference.fa '''\
              '<(path/to/grabix grab file_R1.fastq.gz 1 10000) <(path/to/grabix grab file_R2.fastq.gz 1 10000) | '\
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
        e.assert_any_call('command_bwa', 'command_bwa', 'command_bwa', 'command_bwa', 'command_bwa',
                          cpu=1, job_name='bwa_split', memory=2, working_dir='tests/assets/jobs/test_dataset')
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


class TestMergeBamAndDup(TestVariantCalling):
    """
    Merge bam file generated in bwa mem and defined in SplitFastqStage.
    """

    def test_all_chunk_bam_list_file(self):
        pass

    def test_merge_command(self):
        pass

    def _run(self):
        pass


class TestPostAlignmentScatter(TestVariantCalling):

    def test_split_genome_files(self):
       pass

    def test_split_genome_in_chunks(self):
        pass

    def test_split_genome_chromosomes(self):
        pass


class TestScatterBaseRecalibrator(TestVariantCalling):

    def test_base_recalibrator_cmd(self):
        pass

    def test_run(self):
        pass


class TestGatherBQSRReport(TestVariantCalling):
    def test_run(self):
        pass


class TestScatterApplyBQSR(TestVariantCalling):

    def test_apply_bqsr_cmd(self, chrom_names):
        pass

    def test_run(self):
        pass


class TestGatherRecalBam(TestVariantCalling):

    def test_run(self):
        pass


class TestSplitHaplotypeCaller(TestVariantCalling):

    def test_haplotype_caller_cmd(self, chunk, region_file):
        pass

    def test_run(self):
        pass


class TestGatherGVCF(TestVariantCalling):
    def test_run(self):
        pass


class TestSplitGenotypeGVCFs(TestVariantCalling):

    def test_genotypegvcf_cmd(self):
        pass

    def test_run(self):
        pass


class TestGatherVCF(TestVariantCalling):

    def test_run(self):
        pass


class TestVariantAnnotation(TestVariantCalling):

    def test_run(self):
        pass


class TestVQSRFiltrationSNPs(TestVariantCalling):
    def test_run(self):
        pass


class TestVQSRFiltrationIndels(TestVariantCalling):
    def test_run(self):
        pass


class TestApplyVQSRSNPs(TestVariantCalling):
    def test_run(self):
        pass


class TestApplyVQSRIndels(TestVariantCalling):
    def test_run(self):
        pass


class TestMergeVariants(TestVariantCalling):

    def test_run(self):
        pass


class TestSelectVariants(TestVariantCalling):

    def test_run(self):
        pass

class TestSNPsFiltration():
    def test_run(self):
        pass

class TestIndelsFiltration(TestVariantCalling):
    def test_run(self):
        pass


class TestSplitHaplotypeCaller(TestVariantCalling):

    def setUp(self):
        self.stage = SplitHaplotypeCaller(dataset=self.dataset)

    def test_split_genome_in_chunks(self):
        self.dataset.reference_genome = join(self.assets_path, 'genome.fa')
        print('\n'.join(str(e) for e in self.stage.split_genome_in_chunks()))



