__author__ = 'mwham'
import os
import shutil
from tests.test_analysisdriver import TestAnalysisDriver
from analysis_driver import process_trigger
from analysis_driver.dataset_scanner import RunScanner
from analysis_driver.config import default as cfg


class TestProcessTrigger(TestAnalysisDriver):
    dataset = 'test_dataset'

    @property
    def from_dir(self):
        return os.path.join(self.data_transfer, 'from')

    @property
    def to_dir(self):
        return os.path.join(self.data_transfer, 'to')

    def setUp(self):
        cfg = {
            'lock_file_dir' : os.path.join(self.data_transfer, 'from'),
            'input_dir' : os.path.join(self.data_transfer, 'from')
        }
        self.scanner = RunScanner(cfg)

        self.scanner.reset(self.dataset)

    def tearDown(self):
        for f in os.listdir(os.path.join(self.to_dir, self.dataset)):
            os.remove(os.path.join(self.to_dir, self.dataset, f))


class TestTransfer(TestProcessTrigger):
    def test_transfer(self):
        print(os.getcwd())
        print(os.path.exists("tests/assets/data_transfer/from/"))
        process_trigger._transfer_to_int_dir(self.dataset, self.from_dir, self.to_dir, 1)

        new_dataset = os.path.join(self.to_dir, self.dataset)
        observed = os.listdir(new_dataset)
        expected = ['RTAComplete.txt', 'thang', 'thing']
        print(observed, '\n', expected)
        assert observed == expected

