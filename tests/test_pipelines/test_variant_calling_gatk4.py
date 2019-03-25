import os
from os.path import join, dirname
from unittest.mock import patch

from egcg_core.constants import ELEMENT_LANE, ELEMENT_NB_READS_CLEANED, ELEMENT_RUN_NAME, ELEMENT_RUN_ELEMENT_ID

from analysis_driver.pipelines.variant_calling_gatk4 import SplitBWA

from tests.test_pipelines.test_variant_calling import TestVariantCalling

patch_executor = patch('analysis_driver.pipelines.variant_calling_gatk4.executor.execute')


class TestSplitBWA(TestVariantCalling):

    def setUp(self):
        self.current_wd = os.curdir
        os.chdir(dirname(dirname(__file__)))
        self.stage = SplitBWA(dataset=self.dataset)

    def tearDown(self):
        os.chdir(self.current_wd)

    def test_chunks_from_fastq(self):
        chunks = self.stage.chunks_from_fastq([os.path.join(self.assets_path, 'indexed_fastq_file1.gz')])
        assert chunks == [
            (1, 100000000), (100000001, 200000000), (200000001, 300000000), (300000001, 400000000),
            (400000001, 500000000)
        ]

    def test_bwa_command(self):
        cmd = self.stage.bwa_command(
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

        with patch.object(SplitBWA, '_indexed_fastq_for_run_element', return_value=index_fastq_files), \
             patch_executor as e, patch.object(SplitBWA, 'bwa_command', return_value='command_bwa') as mcommand:
            self.stage._run()
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
