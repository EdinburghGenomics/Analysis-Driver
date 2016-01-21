__author__ = 'mwham'
import os
from unittest.mock import patch
from tests.test_analysisdriver import TestAnalysisDriver
from analysis_driver import transfer_data, executor
from analysis_driver.dataset_scanner import RunDataset
from analysis_driver.writer.bash_commands import rsync_from_to
from tests.fake_rest_api import fake_request, DB, endpoints


class TestProcessTrigger(TestAnalysisDriver):

    @property
    def from_dir(self):
        return os.path.join(self.data_transfer, 'from')

    @property
    def to_dir(self):
        return os.path.join(self.data_transfer, 'to')

    @patch('requests.request', new=fake_request)
    def setUp(self):
        self.setup_db(DB, endpoints)
        self.dataset = RunDataset(
            name='test_dataset',
            path=os.path.join(self.from_dir, 'test_dataset'),
            use_int_dir=False
        )

    def tearDown(self):
        for f in os.listdir(os.path.join(self.to_dir, self.dataset.name)):
            os.remove(os.path.join(self.to_dir, self.dataset.name, f))

    def test_transfer(self):
        transfer_data._transfer_run_to_int_dir(
            self.dataset,
            self.from_dir,
            self.to_dir,
            1,
            rsync_append_verify=False
        )
        new_dataset = os.path.join(self.to_dir, self.dataset.name)
        self.compare_lists(
            observed=os.listdir(new_dataset),
            expected=['RTAComplete.txt', 'thang', 'thing']
        )

    def test_exclude(self):
        cmd = rsync_from_to(
            os.path.join(self.from_dir, self.dataset.name),
            self.to_dir,
            append_verify=False,
            exclude='thing'
        )
        executor.execute([cmd], env='local').join()
        self.compare_lists(
            observed=os.listdir(os.path.join(self.to_dir, self.dataset.name)),
            expected=['RTAComplete.txt', 'thang']
        )
