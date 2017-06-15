from unittest.mock import Mock, patch
from analysis_driver.pipelines.demultiplexing import FastqFilter
from tests.test_analysisdriver import TestAnalysisDriver, NamedMock


class TestFastqFilter(TestAnalysisDriver):
    def test_run(self):
        run_info = Mock(reads=Mock(
            upstream_read=Mock(attrib={'NumCycles': '151'}),
            downstream_read=Mock(attrib={'NumCycles': '151'}),
            index_lengths=[8]
        ))
        dataset = NamedMock(
            real_name='test',
            lane_metrics=[{'pc_q30': 73, 'lane_number': 3}, {'pc_q30': 73, 'lane_number': 4}],
            run_info=run_info
        )
        f = FastqFilter(dataset=dataset)

        results_fastq = [
            [('fastq_L1_R1_001.fastq.gz', 'fastq_L1_R2_001.fastq.gz')],
            [('fastq_L2_R1_001.fastq.gz', 'fastq_L2_R2_001.fastq.gz')],
            [('fastq_L3_R1_001.fastq.gz', 'fastq_L3_R2_001.fastq.gz')],
            [('fastq_L4_R1_001.fastq.gz', 'fastq_L4_R2_001.fastq.gz')],
            [('fastq_L5_R1_001.fastq.gz', 'fastq_L5_R2_001.fastq.gz')],
            [('fastq_L6_R1_001.fastq.gz', 'fastq_L6_R2_001.fastq.gz')],
            [('fastq_L7_R1_001.fastq.gz', 'fastq_L7_R2_001.fastq.gz')],
            [('fastq_L8_R1_001.fastq.gz', 'fastq_L8_R2_001.fastq.gz')]
        ]
        patch_executor = patch('analysis_driver.pipelines.demultiplexing.executor.execute')
        patch_detector = patch('analysis_driver.pipelines.demultiplexing.BadTileCycleDetector')
        patch_find = patch('analysis_driver.pipelines.demultiplexing.find_all_fastq_pairs_for_lane',
                           side_effect=results_fastq + results_fastq)

        with patch_find, patch_executor as pexecute, patch_detector as pdetector:
            instance = pdetector.return_value
            instance.detect_bad_tiles.return_value = {3: [1101]}
            instance.detect_bad_cycles.return_value = {4: [310, 308, 307, 309]}
            f._run()

            expected_call_l3 = (
                'run_filterer fastq_L3_R1_001.fastq.gz fastq_L3_R2_001.fastq.gz '
                'fastq_L3_R1_001_filtered.fastq.gz fastq_L3_R2_001_filtered.fastq.gz '
                'fastq_L3_R1_001_filtered.fastq fastq_L3_R2_001_filtered.fastq '
                'fastq_L3.log --stats_file fastq_L3_fastqfilterer.stats --remove_tiles 1101'
            )
            expected_call_l4 = (
                'run_filterer fastq_L4_R1_001.fastq.gz fastq_L4_R2_001.fastq.gz '
                'fastq_L4_R1_001_filtered.fastq.gz fastq_L4_R2_001_filtered.fastq.gz '
                'fastq_L4_R1_001_filtered.fastq fastq_L4_R2_001_filtered.fastq '
                'fastq_L4.log --stats_file fastq_L4_fastqfilterer.stats --trim_r2 147'
            )
            assert expected_call_l3 == pexecute.call_args[0][2]
            assert expected_call_l4 == pexecute.call_args[0][3]
