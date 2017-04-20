import csv
import time
import illuminate
from os.path import join, isfile
from bitstring import ReadError
from egcg_core import executor
from egcg_core.util import str_join
from analysis_driver.exceptions import AnalysisDriverError
from analysis_driver.reader.run_info import Reads
from analysis_driver.segmentation import Stage


class BCLValidator(Stage):
    @property
    def run_dir(self):
        return self.dataset.input_dir

    @property
    def basecalls_dir(self):
        return join(self.run_dir, 'Data', 'Intensities', 'BaseCalls')

    @property
    def validation_log(self):
        return join(self.run_dir, 'checked_bcls.csv')

    @property
    def validate_expr(self):
        return str_join(
            'function check_bcl { gzip -t ', self.basecalls_dir, '/${1}; x=$?; echo "$1,$x" >> ', self.validation_log, '; }'
        )

    def call_bcl_check(self):
        bcls = self.get_bcls_to_check()
        if bcls:
            self.run_bcl_check(bcls)

    def check_bcls(self):
        while self.dataset.is_sequencing():
            self.call_bcl_check()
            time.sleep(1200)
        # call bcl check again in case the run is finished but not all bcl files have been checked
        self.call_bcl_check()

    def get_bcls_to_check(self):
        """
        Get all present cycles from the run's InterOp, walk backwards through them to find the last completed
        cycle (i.e. has a full set of tiles), and find all bcls for those cycles/tiles that haven't yet been
        checked.
        """
        all_cycles = self._all_cycles_from_interop()
        ncycles = sum(Reads.num_cycles(r) for r in self.dataset.run_info.reads.reads)
        if all_cycles and all_cycles[-1] > ncycles:
            raise AnalysisDriverError(
                'Number of cycles (%s) disagrees with RunInfo (%s)' % (all_cycles[-1], ncycles)
            )

        tile_ids = self.dataset.run_info.tiles

        last_completed_cycle = 0
        for c in sorted(set(all_cycles), reverse=True):
            if all_cycles.count(c) >= len(tile_ids):
                last_completed_cycle = c + 1  # compensate for zero-indexing
                break

        if last_completed_cycle == 0:  # no cycles are complete, so do nothing
            return []

        validated_bcls = []
        if isfile(self.validation_log):
            with open(self.validation_log, 'r') as f:
                validated_bcls = [l.rstrip('\n') for l in f.readlines()]

        bcls_to_check = []
        for c in range(1, last_completed_cycle):
            cycle_id = 'C%s.1' % c
            for t in tile_ids:
                lane = t[0]  # 1_1101
                lane_id = 'L00' + lane
                bcl = join(lane_id, cycle_id, 's_' + t + '.bcl.gz')
                if bcl not in validated_bcls:
                    bcls_to_check.append(bcl)

        self.info('Will check %s bcls up to cycle %s', len(bcls_to_check), last_completed_cycle - 1)
        return bcls_to_check

    def run_bcl_check(self, bcls, slice_size=50, max_job_number=500):
        """
        Run bcl checks through executor.execute. Commands will be collapsed to 50 sequential commands per
        array job so we don't spam the resource manager with 200,000 commands at once.
        """
        max_nb_bcl = max_job_number * slice_size
        for i in range(0, len(bcls), max_nb_bcl):
            tmp_bcls = bcls[i:i + max_nb_bcl]
            sliced_job_array = [
                ['check_bcl %s' % f for f in tmp_bcls[start:start+slice_size]]
                for start in range(0, len(tmp_bcls), slice_size)
            ]

            executor.execute(
                *['\n'.join(cmd_slice) for cmd_slice in sliced_job_array],
                prelim_cmds=[self.validate_expr],
                job_name='bcl_validation',
                working_dir=self.working_dir,
                log_commands=False,
                cpus=1,
                mem=6
            ).join()
        self.info('Finished validation. Check validation log for exit statuses per file.')

    def run_bcl_check_local(self, bcls, parallel=True):
        """Run bcl checks locally through ArrayExecutor."""
        e = executor.ArrayExecutor(['gzip -t ' + join(self.basecalls_dir, f) for f in bcls], stream=parallel)
        e.start()
        e.join()

        with open(self.validation_log, 'a') as f:
            # cmds and e.executors are in the same order, so zip produces corresponding bcls/exit statuses
            for bcl, exit_status in zip(bcls, e.exit_statuses):
                f.write('%s,%s\n' % (join(self.basecalls_dir, bcl), exit_status))

    def read_invalid_files(self):
        with open(self.validation_log, 'r', newline='') as f:
            reader = csv.reader(f, delimiter=',')
            return [bcl for bcl, exit_status in reader if int(exit_status) != 0]

    def _all_cycles_from_interop(self):
        try:
            return illuminate.InteropDataset(self.run_dir).ExtractionMetrics().data['cycle']
        except (illuminate.InteropFileNotFoundError, ReadError):
            self.warning('Cannot load Interop from %s' % self.run_dir)
            return []
