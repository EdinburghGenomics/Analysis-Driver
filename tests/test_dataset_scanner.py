import pytest

__author__ = 'mwham'
__author__ = 'tcezard'
import os
from tests.test_analysisdriver import TestAnalysisDriver
from analysis_driver.dataset_scanner import RunScanner, DATASET_NEW, DATASET_READY, DATASET_PROCESSING, \
    DATASET_PROCESSED_FAIL, DATASET_PROCESSED_SUCCESS, DATASET_ABORTED, RunDataset


class TestRunDataset(TestAnalysisDriver):
    def setUp(self):
        lock_file_dir = os.path.join(self.data_transfer, 'from')
        self.dataset_ready = RunDataset(
            name='that',
            path=os.path.join(self.data_transfer, 'from', 'that'),
            lock_file_dir=lock_file_dir
        )
        self.dataset_not_ready = RunDataset(
            name='this',
            path=os.path.join(self.data_transfer, 'from', 'this'),
            lock_file_dir=lock_file_dir
        )
        self.dataset_ready.reset()
        self.dataset_not_ready.reset()

    def tearDown(self):
        self.dataset_ready.reset()
        self.dataset_not_ready.reset()

    def test_change_status_not_ready(self):
        #dataset not ready
        self.dataset_not_ready.abort()
        assert self.dataset_not_ready.dataset_status == DATASET_ABORTED
        self.dataset_not_ready.reset()
        assert self.dataset_not_ready.dataset_status == DATASET_NEW
        with pytest.raises(AssertionError):
            self.dataset_not_ready.start()
        with pytest.raises(AssertionError):
            self.dataset_not_ready.fail()
        with pytest.raises(AssertionError):
            self.dataset_not_ready.succeed()

    def test_change_status_ready(self):
        self.dataset_ready.reset()
        assert self.dataset_ready.dataset_status == DATASET_READY
        self.dataset_ready.start()
        assert self.dataset_ready.dataset_status == DATASET_PROCESSING
        self.dataset_ready.fail()
        assert self.dataset_ready.dataset_status == DATASET_PROCESSED_FAIL
        self.dataset_ready.reset()
        self.dataset_ready.start()
        self.dataset_ready.succeed()
        assert self.dataset_ready.dataset_status == DATASET_PROCESSED_SUCCESS



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
        self.scanner.get('this').reset()
        self.scanner.get('other').reset()
        more = self.scanner.get('more')
        more.reset()
        more.start()
        that = self.scanner.get('that')
        that.reset()
        that.start()
        that.fail()
        with open(self.triggerignore, 'w') as f:
            for d in ['test_dataset\n']:
                f.write(d)

    def tearDown(self):
        more = self.scanner.get('more')
        more.reset()
        that = self.scanner.get('that')
        that.reset()
        os.remove(self.triggerignore)

    def test_scan_datasets(self):
        datasets = self.scanner.scan_datasets()
        print(datasets)
        print(os.listdir(os.path.join(self.data_transfer, 'from')))
        for d in os.listdir(os.path.join(self.data_transfer, 'from')):
            print(os.listdir(os.path.join(self.data_transfer, 'from', d)))
        for observed, expected in (
            (set([str(s) for s in datasets[DATASET_NEW]]), set(['this', 'other'])),
            (set([str(s) for s in datasets[DATASET_READY]]), set(['another'])),
            (set([str(s) for s in datasets[DATASET_PROCESSING]]), set(['more'])),
            (set([str(s) for s in datasets[DATASET_PROCESSED_FAIL]]), set(['that'])),
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



    @staticmethod
    def _flatten(d):
        l = list()
        for k, v in d.items():
            if type(v) is list:
                l.extend(v)
            else:
                l.append(v)
        return l


