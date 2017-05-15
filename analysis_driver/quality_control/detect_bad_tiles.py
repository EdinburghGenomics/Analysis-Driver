from collections import defaultdict
from itertools import islice

from analysis_driver.reader.run_info import Reads
from .quality_control_base import QualityControl
from interop.py_interop_run_metrics import run_metrics as RunMetrics



class BadTileDetector(QualityControl):
    def __init__(self, job_dir, dataset):
        super().__init__(dataset, job_dir)
        self.run_dir = dataset.input_dir
        self.tile_ids = dataset.run_info.tiles
        self.ncycles = sum(Reads.num_cycles(r) for r in dataset.run_info.reads.reads)
        self.window_size = 50
        self.quality_threshold = 20

        self.read_interop_metrics()

    def read_interop_metrics(self):
        run_metrics = RunMetrics()
        run_metrics.read(self.run_dir)
        q_metric_set = run_metrics.q_metric_set()
        metrics = q_metric_set.metrics()

        self.all_lanes = {}
        for lane in range(1, 9):
            self.all_lanes[lane] = {}

        for i in range(metrics.size()):
            tile_dict = self.all_lanes[metrics[i].lane()]
            if metrics[i].tile() not in tile_dict:
                tile_dict[metrics[i].tile()] = [None] * self.ncycles
            tile_dict[metrics[i].tile()][metrics[i].cycle()-1] = metrics[i].qscore_hist()

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
        d = q8 * 8 + q12 * 12 + q22 * 22 + q27 * 27 + q32 * 32 + q37 * 37 + q41 * 41
        return d/sum((q8, q12, q22, q27, q32, q37, q41))

    def is_bad_sliding_window(self, lane, tile):
        cycles_list = self.all_lanes[lane][tile]
        for window in self.windows(cycles_list, self.window_size):
            avg = self.average_from_list_hist(window)
            if avg is not None and avg < self.quality_threshold:
                return True
        return False

    def detect_bad_tile(self):
        bad_tiles_per_lanes = defaultdict(list)
        for lane in self.all_lanes:
            for tile in self.all_lanes[lane]:
                if self.is_bad_sliding_window(lane, tile):
                    bad_tiles_per_lanes[lane].append(tile)
        return bad_tiles_per_lanes
