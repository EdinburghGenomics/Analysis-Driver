from unittest.mock import Mock, patch

from analysis_driver.dataset import NoCommunicationDataset
from analysis_driver.pipelines.demultiplexing import FastqFilter
from analysis_driver.quality_control.detect_bad_cycles_tiles import BadTileCycleDetector
from tests.test_analysisdriver import TestAnalysisDriver, NamedMock


class TestFastqFilter(TestAnalysisDriver):

    def setUp(self):
        self.dataset = NoCommunicationDataset('test')

    def test_run(self):
        run_info = Mock(reads=Mock(
            upstream_read=Mock(attrib={'NumCycles': '151'}),
            downstream_read=Mock(attrib={'NumCycles': '151'}),
            index_lengths=[8]
        ))
        dataset = NamedMock(
            real_name='test',
            lane_metrics=[{'pc_q30':73, 'lane_number': 3}],
            run_info=run_info
        )
        f=FastqFilter(dataset=dataset)
        results_fastq = [
            [('fastq_L1_R1', 'fastq_L1_R2')],
            [('fastq_L2_R1', 'fastq_L2_R2')],
            [('fastq_L3_R1', 'fastq_L3_R2')],
            [('fastq_L4_R1', 'fastq_L4_R2')],
            [('fastq_L5_R1', 'fastq_L5_R2')],
            [('fastq_L6_R1', 'fastq_L6_R2')],
            [('fastq_L7_R1', 'fastq_L7_R2')],
            [('fastq_L8_R1', 'fastq_L8_R2')]
        ]
        patch_find = patch('analysis_driver.pipelines.demultiplexing.find_all_fastq_pairs_for_lane',
                           side_effect=results_fastq)
        patch_executor = patch('analysis_driver.pipelines.demultiplexing.executor.execute')
        patch_detector = patch('analysis_driver.pipelines.demultiplexing.BadTileCycleDetector')
        with patch_find, patch_executor as pexecute, patch_detector as pdetector:
            instance = pdetector.return_value
            instance.detect_bad_tile.return_value = {3: ['1101']}
            f._run()
            expected_call_L3 = """mkfifo fastq_L3_R1_filtered.fastq
mkfifo fastq_L3_R2_filtered.fastq
set -e; path/to/fastq-filterer --i1 fastq_L3_R1 --i2 fastq_L3_R2 --o1 fastq_L3_R1_filtered.fastq --o2 fastq_L3_R2_filtered.fastq --threshold 36 --remove_tiles 1101 & pigz -c -p 10 fastq_L3_R1_filtered.fastq > fastq_L3_R1_filtered.fastq.gz & pigz -c -p 10 fastq_L3_R2_filtered.fastq > fastq_L3_R2_filtered.fastq.gz
EXIT_CODE=$?
rm fastq_L3_R1_filtered.fastq fastq_L3_R2_filtered.fastq
(exit $EXIT_CODE) && mv fastq_L3_R1_filtered.fastq.gz fastq_L3_R1
(exit $EXIT_CODE) && mv fastq_L3_R2_filtered.fastq.gz fastq_L3_R2
(exit $EXIT_CODE)"""
            assert expected_call_L3 == pexecute.call_args[0][2]
