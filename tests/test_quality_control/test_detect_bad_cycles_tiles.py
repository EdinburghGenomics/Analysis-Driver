import os
from unittest.mock import Mock, patch, call

from analysis_driver.quality_control.detect_bad_cycles_tiles import BadTileCycleDetector
from tests.test_analysisdriver import TestAnalysisDriver


class TestBadTileCycleDetector(TestAnalysisDriver):
    def setUp(self):
        run_info = Mock(
            tiles=('1_1101', '2_1101', '1_1102', '2_1102'),
            reads=Mock(reads=[Mock(attrib={'NumCycles': '3'})])
        )
        self.job_dir = os.path.join(TestAnalysisDriver.assets_path, 'bcl_validation')
        self.detector = BadTileCycleDetector(Mock(input_dir=self.job_dir, run_info=run_info))

    def test_windows(self):
        l = ['a', 'b', 'c', 'd', 'e', 'f', 'g']
        expected_lists = [('a', 'b', 'c'), ('b', 'c', 'd'), ('c', 'd', 'e'), ('d', 'e', 'f'), ('e', 'f', 'g')]
        assert list(self.detector.windows(l, window_size=3)) == expected_lists
        l = ['a', 'b', 'c']
        assert list(self.detector.windows(l, window_size=4)) == []


    def test_average_from_list_hist(self):
        list_hist = [
            (1, 2, 1, 5, 6, 9, 3),
            (0, 1, 2, 4, 1, 4, 8),
        ]
        assert self.detector.average_from_list_hist(list_hist) == 32.1063829787234
        list_hist = [(1, 1, 1, 1, 1, 1, 1)]
        assert self.detector.average_from_list_hist(list_hist) == 25.571428571428573
        list_hist = [(1, 1, 1, 1, 1, 1, 1), None]
        assert self.detector.average_from_list_hist(list_hist) == 25.571428571428573
        list_hist = [None, None]
        assert self.detector.average_from_list_hist(list_hist) is None
        list_hist = []
        assert self.detector.average_from_list_hist(list_hist) is None

    def test_is_bad_tile_sliding_window(self):
        list_q_hist1 = [
            (0, 1, 1, 1, 1, 1, 10),
            (0, 1, 1, 1, 1, 1, 10),
            (0, 10, 1, 1, 1, 1, 1),
            (0, 10, 1, 1, 1, 1, 1),
            (0, 1, 1, 1, 1, 1, 10)
        ]
        list_q_hist2 = [
            (0, 1, 1, 1, 1, 1, 10),
            (0, 1, 1, 1, 1, 1, 10),
            (0, 10, 1, 1, 1, 1, 1),
            (0, 1, 1, 1, 1, 1, 10),
            (0, 1, 1, 1, 1, 1, 10)
        ]
        self.detector.window_size = 2

        assert self.detector.is_bad_tile_sliding_window(list_q_hist1)

        assert not self.detector.is_bad_tile_sliding_window(list_q_hist2)

    def test_detect_bad_tile(self):
        list_q_hist1 = [
            (0, 1, 1, 1, 1, 1, 10),
            (0, 1, 1, 1, 1, 1, 10),
            (0, 10, 1, 1, 1, 1, 1),
            (0, 10, 1, 1, 1, 1, 1),
            (0, 1, 1, 1, 1, 1, 10)
        ]
        list_q_hist2 = [
            (0, 1, 1, 1, 1, 1, 10),
            (0, 1, 1, 1, 1, 1, 10),
            (0, 10, 1, 1, 1, 1, 1),
            (0, 1, 1, 1, 1, 1, 10),
            (0, 1, 1, 1, 1, 1, 10)
        ]
        self.detector.window_size = 2
        self.detector.all_lanes = {1: ({'1101': list_q_hist1, '1102': list_q_hist2}, {})}
        assert dict(self.detector.detect_bad_tile()) == {1: ['1101']}

    def test_detect_bad_cycles(self):
        list_q_hist1 = [
            (0, 10, 1, 1, 1, 1, 1),
            (0, 10, 1, 1, 1, 1, 1),
            (0, 10, 1, 1, 1, 1, 1),
            (0, 10, 1, 1, 1, 1, 1),
            (10, 1, 1, 1, 1, 1, 1)
        ]
        list_q_hist2 = [
            (0, 1, 1, 1, 1, 1, 10),
            (0, 1, 1, 1, 1, 1, 10),
            (0, 10, 1, 1, 1, 1, 1),
            (0, 1, 1, 1, 1, 1, 10),
            (0, 1, 1, 1, 1, 1, 10)
        ]
        self.detector.all_lanes = {1: ({}, {'1': list_q_hist1, '2': list_q_hist2})}
        assert dict(self.detector.detect_bad_cycle()) == {1: ['1']}
