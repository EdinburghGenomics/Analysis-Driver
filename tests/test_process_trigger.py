__author__ = 'mwham'
import os
import shutil
from tests.test_analysisdriver import TestAnalysisDriver
from analysis_driver import process_trigger
from analysis_driver.dataset_scanner import RunScanner
from analysis_driver.config import default as cfg


class TestProcessTrigger(TestAnalysisDriver):
    dataset = 'test_dataset'
    scanner = RunScanner(cfg)
    @property
    def from_dir(self):
        return os.path.join(self.data_transfer, 'from')

    @property
    def to_dir(self):
        return os.path.join(self.data_transfer, 'to')

    @property
    def triggerignore(self):
        return os.path.join(self.from_dir, '.triggerignore')

    def setUp(self):
        self.scanner.reset(self.dataset)

    def tearDown(self):
        with open(os.path.join(self.from_dir, '.triggerignore'), 'w') as f:
            for d in ['dir_t?_be_ign*d\n', 'test_dataset\n']:
                f.write(d)


class TestTransfer(TestProcessTrigger):
    def test_transfer(self):
        self.scanner.reset(self.dataset)
        print(os.getcwd())
        print(os.path.exists("tests/assets/data_transfer/from/"))
        process_trigger._transfer_to_int_dir(self.dataset, self.from_dir, self.to_dir, 1)

        new_dataset = os.path.join(self.to_dir, self.dataset)
        observed = os.listdir(new_dataset)
        expected = ['RTAComplete.txt', 'thang', 'thing']
        print(observed, '\n', expected)
        assert observed == expected

