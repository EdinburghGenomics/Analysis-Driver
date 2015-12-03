__author__ = 'mwham'
import analysis_driver
import os
import shutil
from tests.test_analysisdriver import TestAnalysisDriver
from analysis_driver import transfer_data
from analysis_driver.dataset_scanner import RunDataset


class TestProcessTrigger(TestAnalysisDriver):

    @property
    def from_dir(self):
        return os.path.join(self.data_transfer, 'from')

    @property
    def to_dir(self):
        return os.path.join(self.data_transfer, 'to')

    def setUp(self):
        self.dataset = RunDataset(
            name='test_dataset',
            path=os.path.join(self.from_dir, 'test_dataset'),
            lock_file_dir=self.from_dir
        )
        self.dataset.reset()

    def tearDown(self):
        for f in os.listdir(os.path.join(self.to_dir, self.dataset.name)):
            os.remove(os.path.join(self.to_dir, self.dataset.name, f))
        self.dataset.reset()

    def test_transfer(self):
        transfer_data._transfer_to_int_dir(
            self.dataset,
            self.from_dir,
            self.to_dir,
            1,
            rsync_append_verify=False
        )
        new_dataset = os.path.join(self.to_dir, self.dataset.name)
        observed = os.listdir(new_dataset)
        expected = ['RTAComplete.txt', 'thang', 'thing']
        print(observed, '\n', expected)
        assert observed == expected
