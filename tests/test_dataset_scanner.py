import pytest

__author__ = 'mwham'
__author__ = 'tcezard'
import os
from tests.test_analysisdriver import TestAnalysisDriver
from analysis_driver.dataset_scanner import RunScanner, DATASET_NEW, DATASET_READY, DATASET_PROCESSING, \
    DATASET_PROCESSED_FAIL, DATASET_PROCESSED_SUCCESS, DATASET_ABORTED

class TestRunScanner(TestAnalysisDriver):

    @property
    def triggerignore(self):
        return os.path.join(self.scanner.lock_file_dir, '.triggerignore')

    def setUp(self):
        cfg = {
            'lock_file_dir' : os.path.join(self.data_transfer, 'from'),
            'input_dir' : os.path.join(self.data_transfer, 'from')
        }
        self.scanner = RunScanner(cfg)
        self.scanner.reset('this')
        self.scanner.reset('that')
        self.scanner.start('that')
        self.scanner.fail('that')
        self.scanner.reset('more')
        self.scanner.start('more')
        with open(self.triggerignore, 'w') as f:
            for d in ['test_dataset\n']:
                f.write(d)

    def tearDown(self):
        self.scanner.reset('more')
        self.scanner.reset('that')
        os.remove(self.triggerignore)

    def test_scan_datasets(self):
        datasets = self.scanner.scan_datasets()
        print(datasets)

        for observed, expected in (
            (set(datasets[DATASET_NEW]), set(['this', 'other'])),
            (set(datasets[DATASET_READY]), set(['another'])),
            (set(datasets[DATASET_PROCESSING]), set(['more'])),
            (set(datasets[DATASET_PROCESSED_FAIL]), set(['that'])),
        ):
            assert observed == expected

    def test_triggerignore(self):
        with open(self.triggerignore, 'r') as f:
            assert f.readlines() == ['test_dataset\n']

        expected = ['other\n', 'that\n', 'this\n']
        with open(self.triggerignore, 'w') as f:
            for d in expected:
                f.write(d)

        with open(self.triggerignore, 'r') as f:
            assert f.readlines() == expected
            for d in expected:
                assert d.strip() not in self._flatten(self.scanner.scan_datasets())

    def test_switch_status(self):
        #dataset not ready
        d = 'this'
        self.scanner.abort(d)
        assert self.scanner.dataset_status(d) == DATASET_ABORTED
        self.scanner.reset(d)
        assert self.scanner.dataset_status(d) == DATASET_NEW
        with pytest.raises(AssertionError):
            self.scanner.start(d)
            self.scanner.fail(d)
            self.scanner.succeed(d)

        d = 'that'
        self.scanner.reset(d)
        assert self.scanner.dataset_status(d) == DATASET_READY
        self.scanner.start(d)
        assert self.scanner.dataset_status(d) == DATASET_PROCESSING
        self.scanner.fail(d)
        assert self.scanner.dataset_status(d) == DATASET_PROCESSED_FAIL
        self.scanner.reset(d)
        self.scanner.start(d)
        self.scanner.succeed(d)
        assert self.scanner.dataset_status(d) == DATASET_PROCESSED_SUCCESS


    @staticmethod
    def _flatten(d):
        l = list()
        for k, v in d.items():
            if type(v) is list:
                l.extend(v)
            else:
                l.append(v)
        return l


