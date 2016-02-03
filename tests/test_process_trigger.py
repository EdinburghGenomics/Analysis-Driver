__author__ = 'mwham'
import os
import shutil
from unittest.mock import patch
from tests.test_analysisdriver import TestAnalysisDriver
from analysis_driver import transfer_data
from analysis_driver.dataset_scanner import RunDataset
from analysis_driver.util.bash_commands import rsync_from_to
from analysis_driver.config import default as cfg
from tests.test_dataset_scanner import patched_request


class TestProcessTrigger(TestAnalysisDriver):

    @property
    def from_dir(self):
        return os.path.join(self.data_transfer, 'from')

    @property
    def to_dir(self):
        return os.path.join(self.data_transfer, 'to')

    def _patched_rsync(self, from_dir, to_dir, exclude=None):
        cmd = rsync_from_to(from_dir + '/', to_dir, exclude=exclude)
        cmd = cmd.replace('--append-verify ', '')
        return patch('analysis_driver.transfer_data.rsync_from_to', return_value=cmd)

    def setUp(self):
        with patched_request:
            self.dataset = RunDataset(
                name='test_dataset',
                path=os.path.join(self.from_dir, 'test_dataset'),
                use_int_dir=False
            )
        self.old_job_dir = cfg['jobs_dir']
        cfg.content['jobs_dir'] = self.data_transfer
        os.makedirs(os.path.join(self.data_transfer, 'test_dataset'), exist_ok=True)

    def tearDown(self):
        for f in os.listdir(os.path.join(self.to_dir, self.dataset.name)):
            os.remove(os.path.join(self.to_dir, self.dataset.name, f))
        shutil.rmtree(os.path.join(self.data_transfer, 'test_dataset'))
        cfg.content['jobs_dir'] = self.old_job_dir

    def test_transfer(self):
        with self._patched_rsync(self.from_dir, self.to_dir):
            transfer_data._transfer_run_to_int_dir(
                self.dataset,
                self.from_dir,
                self.to_dir,
                1
            )
        new_dataset = os.path.join(self.to_dir, self.dataset.name)
        self.compare_lists(
            observed=os.listdir(new_dataset),
            expected=['RTAComplete.txt', 'thang', 'thing']
        )

    def test_exclude(self):
        cmd = rsync_from_to(
            os.path.join(self.from_dir, self.dataset.name, ''),
            self.to_dir,
            exclude='thing'
        )
        assert '--exclude=thing' in cmd
