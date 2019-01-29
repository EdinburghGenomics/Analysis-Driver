import os
import csv
import time
from egcg_core import executor
from egcg_core.util import str_join
from analysis_driver.exceptions import AnalysisDriverError
from analysis_driver.quality_control import interop_metrics
from analysis_driver.reader.run_info import Reads
from analysis_driver.segmentation import Stage


class BCLValidator(Stage):
    @property
    def run_dir(self):
        return self.dataset.input_dir

    @property
    def basecalls_dir(self):
        return os.path.join(self.run_dir, 'Data', 'Intensities', 'BaseCalls')

    @property
    def validation_log(self):
        return os.path.join(self.run_dir, 'checked_bcls.csv')

    def validate_expr(self, validation_log):
        return str_join(
            'function check_bcl { gzip -t ',
            self.basecalls_dir,
            '/${1}; x=$?; echo "$1,$x" >> ',
            validation_log,
            '; }'
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
            raise AnalysisDriverError('Number of cycles (%s) disagrees with RunInfo (%s)' % (all_cycles[-1], ncycles))

        tile_ids = self.dataset.run_info.tiles
        last_completed_cycle = 0
        for c in sorted(set(all_cycles), reverse=True):
            if all_cycles.count(c) >= len(tile_ids):
                last_completed_cycle = c + 1  # compensate for zero-indexing
                break

        if last_completed_cycle == 0:  # no cycles are complete, so do nothing
            return []

        validated_bcls = self.read_valid_files()
        bcls_to_check = []
        for c in range(1, last_completed_cycle):
            cycle_id = 'C%s.1' % c
            for t in tile_ids:
                lane = t[0]  # 1_1101
                lane_id = 'L00' + lane
                bcl = os.path.join(lane_id, cycle_id, 's_' + t + '.bcl.gz')
                if bcl not in validated_bcls:
                    bcls_to_check.append(bcl)

        self.info('Will check %s bcls up to cycle %s', len(bcls_to_check), last_completed_cycle - 1)
        return bcls_to_check

    def run_bcl_check(self, bcls, slice_size=100, max_job_number=500):
        """
        Run bcl checks through executor.execute. Commands will be collapsed to 100 sequential commands per
        array job so we don't spam the resource manager with 200,000 commands at once.
        :param list bcls: File paths of bcls to check
        :param int slice_size: Group size for collapsing array jobs
        :param int max_job_number: Maximum number of bcls that will be processed for this method call
        """
        validation_log_tmp = os.path.join(self.job_dir, 'tmp_checked_bcls.csv')
        if os.path.isfile(validation_log_tmp):
            os.remove(validation_log_tmp)

        max_nb_bcl = max_job_number * slice_size
        for i in range(0, len(bcls), max_nb_bcl):
            tmp_bcls = bcls[i:i + max_nb_bcl]
            sliced_job_array = [
                ['check_bcl %s' % f for f in tmp_bcls[start:start+slice_size]]
                for start in range(0, len(tmp_bcls), slice_size)
            ]

            executor.execute(
                *['\n'.join(cmd_slice) for cmd_slice in sliced_job_array],
                prelim_cmds=[self.validate_expr(validation_log_tmp)],
                job_name='bcl_validation',
                working_dir=self.job_dir,
                log_commands=False,
                cpus=1,
                mem=6
            ).join()

        valid_bcls = self.read_valid_files()
        # Merge the valid bcls and the tested bcls
        with open(self.validation_log, 'w') as f:
            for bcl in valid_bcls:
                f.write('%s,0\n' % bcl)
            if os.path.isfile(validation_log_tmp):
                for bcl, exit_status in self.read_check_bcl_files(validation_log_tmp):
                    f.write('%s,%s\n' % (bcl, exit_status))
        self.info('Finished validation. Found %s invalid files' % len(self.read_invalid_files()))

    def read_check_bcl_files(self, validation_log=None):
        with open(validation_log or self.validation_log, 'r', newline='') as f:
            return [i for i in csv.reader(f, delimiter=',')]

    def read_invalid_files(self):
        if os.path.isfile(self.validation_log):
            return [bcl for bcl, exit_status in self.read_check_bcl_files() if int(exit_status) != 0]
        else:
            return []

    def read_valid_files(self):
        if os.path.isfile(self.validation_log):
            return [bcl for bcl, exit_status in self.read_check_bcl_files() if int(exit_status) == 0]
        else:
            return []

    def _all_cycles_from_interop(self):
        return interop_metrics.get_cycles_extracted(self.run_dir)
