__author__ = 'tcezard'
import shutil
import pytest
from unittest.mock import patch
import os
import requests
import subprocess
from tests.test_analysisdriver import TestAnalysisDriver
from tests.fake_rest_api import fake_request, DB, endpoints
from analysis_driver.dataset_scanner import RunScanner, DATASET_NEW, DATASET_READY, DATASET_PROCESSING, \
    DATASET_PROCESSED_FAIL, DATASET_PROCESSED_SUCCESS, DATASET_ABORTED, DATASET_REPROCESS, RunDataset
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
    @patch('requests.request', new=fake_request)
    def setUp(self):
        self.setup_db(DB, endpoints)
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

    @patch('requests.request', new=fake_request)
    def tearDown(self):
        clean(self.base_dir)

    @patch('requests.request', new=fake_request)
    def test_change_status_not_ready(self):
        d = self.dataset_not_ready

        d.abort()
        assert d.dataset_status == DATASET_ABORTED
        d.reset()
        assert d.dataset_status == DATASET_REPROCESS
        d.start()
        assert d.dataset_status == DATASET_PROCESSING
        d.reset()
        with pytest.raises(AssertionError):
            d.fail()

    @patch('requests.request', new=fake_request)
    def test_change_status_ready(self):
        d = self.dataset_ready

        d.reset()
        assert d.dataset_status in (DATASET_READY, DATASET_REPROCESS)
        d.start()
        assert d.dataset_status == DATASET_PROCESSING
        d.fail()
        assert d.dataset_status == DATASET_PROCESSED_FAIL
        d.reset()
        assert d.dataset_status in (DATASET_READY, DATASET_REPROCESS)
        d.start()
        assert d.dataset_status == DATASET_PROCESSING
        d.succeed()
        assert d.dataset_status == DATASET_PROCESSED_SUCCESS


class TestRunScanner(TestAnalysisDriver):
    @property
    def triggerignore(self):
        return os.path.join(self.scanner.lock_file_dir, '.triggerignore')

    @patch('requests.request', new=fake_request)
    def setUp(self):
        self.setup_db(DB, endpoints)
        self.db = DB
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

    @patch('requests.request', new=fake_request)
    def tearDown(self):
        clean(self.base_dir)
        safe_delete('runs', 'analysis_driver_procs')

    @patch('requests.request', new=fake_request)
    def test_scan_datasets(self):
        datasets = self.scanner.scan_datasets()
        print('datasets: ', datasets)

        for status, expected in (
            (DATASET_NEW, {'other', 'dataset_not_ready'}),
            (DATASET_READY, {'dataset_ready', 'another'}),
            (DATASET_PROCESSING, {'more'}),
            (DATASET_PROCESSED_FAIL, {'that'})
        ):
            self.compare_lists(observed=self._dataset_names(datasets, status), expected=expected)

    @patch('requests.request', new=fake_request)
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
