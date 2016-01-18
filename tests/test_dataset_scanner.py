__author__ = 'tcezard'
import shutil
import pytest
import os
import requests
import subprocess
from tests.test_analysisdriver import TestAnalysisDriver
from analysis_driver.dataset_scanner import RunScanner, DATASET_NEW, DATASET_READY, DATASET_PROCESSING, \
    DATASET_PROCESSED_FAIL, DATASET_PROCESSED_SUCCESS, DATASET_ABORTED, RunDataset
from analysis_driver.config import default as cfg


def seed_directories(base_dir):
    directories_to_create = ['another', 'that', 'more', 'other', 'dataset_ready', 'dataset_not_ready']
    file_to_touch = ['another/RTAComplete.txt', 'more/RTAComplete.txt', 'dataset_ready/RTAComplete.txt']
    for d in directories_to_create:
        os.makedirs(os.path.join(base_dir, d), exist_ok=True)
    for f in file_to_touch:
        open(os.path.join(base_dir, f), 'w').close()


def clean(base_dir):
    shutil.rmtree(os.path.join(base_dir))


def mongod_running():
    p = subprocess.Popen(['ps'], stdout=subprocess.PIPE)
    out, err = p.communicate()
    ps = out.decode('utf-8').split('\n')
    for line in ps:
        if 'mongod' in line:
            return True
    return False


def safe_delete(*endpoints):
    if mongod_running():
        for e in endpoints:
            requests.delete(cfg.query('rest_api', 'url').rstrip('/') + '/' + e)
    else:
        print('No local test database is running!')


class TestRunDataset(TestAnalysisDriver):
    def setUp(self):
        safe_delete('runs', 'samples', 'analysis_driver_procs')
        self.base_dir = os.path.join(self.assets_path, 'dataset_scanner')
        seed_directories(self.base_dir)
        self.dataset_ready = RunDataset(
            name='dataset_ready',
            path=os.path.join(self.base_dir, 'dataset_ready'),
            lock_file_dir=self.base_dir,
            use_int_dir=False
        )
        self.dataset_ready.start()

        self.dataset_not_ready = RunDataset(
            name='dataset_not_ready',
            path=os.path.join(self.data_transfer, 'dataset_not_ready'),
            lock_file_dir=self.base_dir,
            use_int_dir=False
        )
        self.dataset_not_ready.start()

    def tearDown(self):
        clean(self.base_dir)
        self.dataset_ready.reset()
        self.dataset_not_ready.reset()

    def test_change_status_not_ready(self):
        d = self.dataset_not_ready

        d.abort()
        assert d.dataset_status == DATASET_ABORTED
        d.reset()
        assert d.dataset_status == DATASET_NEW or d._most_recent_proc().get('rerun')
        d.start()
        assert d.dataset_status == DATASET_PROCESSING
        d.reset()
        with pytest.raises(AssertionError):
            d.fail()
        with pytest.raises(AssertionError):
            d.succeed()

    def test_change_status_ready(self):
        d = self.dataset_ready

        d.reset()
        assert d.dataset_status == DATASET_READY or d._most_recent_proc().get('rerun')
        d.start()
        assert d.dataset_status == DATASET_PROCESSING
        d.fail()
        assert d.dataset_status == DATASET_PROCESSED_FAIL
        d.reset()
        assert d.dataset_status == DATASET_READY or d._most_recent_proc().get('rerun')
        d.start()
        assert d.dataset_status == DATASET_PROCESSING
        d.succeed()
        assert d.dataset_status == DATASET_PROCESSED_SUCCESS


class TestRunScanner(TestAnalysisDriver):
    @property
    def triggerignore(self):
        return os.path.join(self.scanner.lock_file_dir, '.triggerignore')

    def setUp(self):
        self.base_dir = os.path.join(self.assets_path, 'dataset_scanner')
        seed_directories(self.base_dir)

        cfg = {
            'lock_file_dir': self.base_dir,
            'input_dir': self.base_dir
        }
        self.scanner = RunScanner(cfg)
        more = self.scanner.get('more')
        more.start()
        that = self.scanner.get('that')
        that.start()
        that.fail()

        with open(self.triggerignore, 'w') as f:
            for d in ['test_dataset\n']:
                f.write(d)

    def tearDown(self):
        clean(self.base_dir)
        safe_delete('runs', 'analysis_driver_procs')

    def test_scan_datasets(self):
        datasets = self.scanner.scan_datasets()
        print('datasets: ', datasets)
        print(os.listdir(self.base_dir))
        for d in os.listdir(self.base_dir):
            if os.path.isdir(os.path.join(self.base_dir, d)):
                print(os.listdir(os.path.join(self.base_dir, d)))
        for status, expected in (
            (DATASET_NEW, {'other', 'another', 'dataset_not_ready'}),
            (DATASET_READY, {'dataset_ready'}),
            (DATASET_PROCESSING, {'more'}),
            (DATASET_PROCESSED_FAIL, {'that'})
        ):
            self.compare_lists(observed=self._dataset_names(datasets, status), expected=expected)

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

    @staticmethod
    def _dataset_names(datasets, status):
        return set([str(s).split(' ')[0] for s in datasets[status]])
