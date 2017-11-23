from collections import defaultdict
from itertools import islice
from egcg_core.config import cfg
from egcg_core.app_logging import AppLogger
from interop.py_interop_run_metrics import run_metrics as RunMetrics
from analysis_driver.reader.run_info import Reads


class BadTileCycleDetector(AppLogger):
    def __init__(self, dataset, window_size=60, tile_quality_threshold=None, cycle_quality_threshold=None):
        self.dataset = dataset
        self.run_dir = dataset.input_dir
        self.tile_ids = dataset.run_info.tiles
        self.ncycles = sum(Reads.num_cycles(r) for r in dataset.run_info.reads.reads)
        self.window_size = window_size
        self.all_lanes = None
        self.tile_quality_threshold = tile_quality_threshold or cfg.query(
            'fastq_filterer', 'tile_quality_threshold', ret_default=20
        )
        self.cycle_quality_threshold = cycle_quality_threshold or cfg.query(
            'fastq_filterer', 'cycle_quality_threshold', ret_default=18
        )

        self.read_interop_metrics()

    def read_interop_metrics(self):
        run_metrics = RunMetrics()
        run_metrics.read(self.run_dir)
        q_metric_set = run_metrics.q_metric_set()
        metrics = q_metric_set.metrics()

        self.all_lanes = {}
        for lane in range(1, 9):
            self.all_lanes[lane] = ({}, {})

        for i in range(metrics.size()):
            tile_dict, cycle_dict = self.all_lanes[metrics[i].lane()]
            if metrics[i].tile() not in tile_dict:
                tile_dict[metrics[i].tile()] = [None] * self.ncycles
            if metrics[i].cycle() not in cycle_dict:
                cycle_dict[metrics[i].cycle()] = []

            tile_dict[metrics[i].tile()][metrics[i].cycle()-1] = metrics[i].qscore_hist()
            cycle_dict[metrics[i].cycle()].append(metrics[i].qscore_hist())

    @staticmethod
    def windows(l, window_size):
        it = iter(l)
        result = tuple(islice(it, window_size))
        if len(result) == window_size:
            yield result
        for elem in it:
            result = result[1:] + (elem,)
            yield result

    @staticmethod
    def average_from_list_hist(q_hist_list):
        q_hist_list = [hist for hist in q_hist_list if hist]
        if len(q_hist_list) == 0:
            return None
        q8, q12, q22, q27, q32, q37, q41 = q_hist_list[0]
        for q_hist in q_hist_list[1:]:
            tq8, tq12, tq22, tq27, tq32, tq37, tq41 = q_hist
            q8 += tq8
            q12 += tq12
            q22 += tq22
            q27 += tq27
            q32 += tq32
            q37 += tq37
            q41 += tq41
        n = q8 * 8 + q12 * 12 + q22 * 22 + q27 * 27 + q32 * 32 + q37 * 37 + q41 * 41
        d = sum((q8, q12, q22, q27, q32, q37, q41))
        if d > 0:
            return n / d
        return None

    def is_bad_tile_sliding_window(self, lane, tile, cycles_list):
        window_count = 0
        for window in self.windows(cycles_list, self.window_size):
            window_count += 1
            avg = self.average_from_list_hist(window)
            if avg is not None and avg < self.tile_quality_threshold:
                self.info(
                    'lane %s tile %s window %s-%s: average quality %s < %s',
                    lane, tile, window_count, window_count + self.window_size, avg,
                    self.tile_quality_threshold
                )
                return True
        return False

    def detect_bad_tiles(self):
        bad_tiles_per_lanes = defaultdict(list)
        for lane in self.all_lanes:
            tile_dict, cycle_dict = self.all_lanes[lane]
            for tile in tile_dict:
                if self.is_bad_tile_sliding_window(lane, tile, tile_dict[tile]):
                    bad_tiles_per_lanes[lane].append(tile)
        return bad_tiles_per_lanes

    def detect_bad_cycles(self):
        bad_cycle_per_lanes = defaultdict(list)
        for lane in self.all_lanes:
            tile_dict, cycle_dict = self.all_lanes[lane]
            for cycle in cycle_dict:
                avg = self.average_from_list_hist(cycle_dict[cycle])
                if avg < self.cycle_quality_threshold:
                    self.info('Lane %s cycle %s: average quality %s < %s', lane, cycle, avg, self.cycle_quality_threshold)
                    bad_cycle_per_lanes[lane].append(cycle)
        return bad_cycle_per_lanes
