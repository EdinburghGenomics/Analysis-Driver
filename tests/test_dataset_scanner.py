__author__ = 'tcezard'
import shutil
import pytest
from unittest.mock import patch, Mock
import os
import re
from tests.test_analysisdriver import TestAnalysisDriver
from analysis_driver.dataset_scanner import RunScanner, DATASET_NEW, DATASET_READY, DATASET_PROCESSING, \
    DATASET_PROCESSED_FAIL, DATASET_PROCESSED_SUCCESS, DATASET_ABORTED, DATASET_REPROCESS, DATASET_FORCE_READY,\
    RunDataset
from analysis_driver.config import default as cfg


def seed_directories(base_dir):
    directories_to_create = (
        'failed_dataset', 'processing_dataset', 'dataset_ready', 'dataset_not_ready', 'ignored_dataset'
    )
    rta_completes = ('processing_dataset', 'dataset_ready')
    for d in directories_to_create:
        os.makedirs(os.path.join(base_dir, d), exist_ok=True)
        if d in rta_completes:
            open(os.path.join(base_dir, d, 'RTAComplete.txt'), 'w').close()


def clean(base_dir):
    shutil.rmtree(os.path.join(base_dir))


class FakeRestResponse(Mock):
    def json(self):
        return self.content


def _fake_analysis_driver_proc(d_type, d_name, status=None):
    proc = {
        'proc_id': '_'.join((d_type, d_name)),
        'dataset_type': d_type,
        'dataset_name': d_name
    }
    if status:
        proc['status'] = status
    return proc


def fake_scan_request(req_type, url):
    assert req_type == 'GET'
    endpoint, query_string = url.split('/')[-1].split('?')

    dataset_type = re.search(r'"dataset_type":"(.+?)"', query_string).group(1)
    dataset_name = re.search(r'"dataset_name":"(.+?)"', query_string).group(1)

    s = dataset_name.replace('_dataset', '')

    fake_content = _fake_analysis_driver_proc(dataset_type, dataset_name)
    if s in (DATASET_FORCE_READY, DATASET_PROCESSING, DATASET_PROCESSED_SUCCESS, DATASET_PROCESSED_FAIL, DATASET_ABORTED, DATASET_REPROCESS):
        fake_content['status'] = s
    return FakeRestResponse(content={'data': [fake_content]})


empty_proc_data = {'data': []}


def empty_analysis_driver_procs():
    return patch('requests.request', return_value=FakeRestResponse(content={'data': []}))


def fake_analysis_driver_proc(dataset, status):
    return patch(
        'requests.request',
        return_value=FakeRestResponse(
            content={'data': [_fake_analysis_driver_proc(dataset.type, dataset.name, status=status)]}
        )
    )


class TestRunDataset(TestAnalysisDriver):
    def setUp(self):
        self.base_dir = os.path.join(self.assets_path, 'dataset_scanner')
        seed_directories(self.base_dir)

    def tearDown(self):
        clean(self.base_dir)

    def _query_most_recent_proc(self, dataset_name):
        return ''.join(
            (
                cfg.query('rest_api', 'url'),
                '/analysis_driver_procs',
                '?where={"dataset_type":"run","dataset_name":"',
                dataset_name,
                '"}&sort=-_created'
            )
        )

    @empty_analysis_driver_procs()
    def test_dataset_not_ready(self, mocked_instance):
        d = RunDataset(
            name='dataset_not_ready',
            path=os.path.join(self.data_transfer, 'dataset_not_ready'),
            use_int_dir=False
        )
        assert mocked_instance.call_count == 1  # setting a breakpoint here breaks test, makes call_count 2
        assert d.dataset_status == DATASET_NEW
        assert mocked_instance.call_count == 2
        mocked_instance.assert_called_with('GET', self._query_most_recent_proc('dataset_not_ready'))

    @empty_analysis_driver_procs()
    def test_dataset_ready(self, mocked_instance):
        d = RunDataset(
            name='dataset_ready',
            path=os.path.join(self.base_dir, 'dataset_ready'),
            use_int_dir=False
        )
        assert mocked_instance.call_count == 1
        assert d.dataset_status == DATASET_READY
        assert mocked_instance.call_count == 2
        mocked_instance.assert_called_with('GET', self._query_most_recent_proc('dataset_ready'))

    def test_dataset_status(self):
        with empty_analysis_driver_procs():
            d = RunDataset(
                name='test',
                path=os.path.join(self.base_dir, 'test'),
                use_int_dir=False
            )
            assert d.dataset_status == DATASET_NEW
        for expected_status in (
            DATASET_NEW,
            DATASET_PROCESSING,
            DATASET_REPROCESS,
            DATASET_PROCESSED_FAIL,
            DATASET_PROCESSED_SUCCESS
        ):
            with fake_analysis_driver_proc(d, expected_status):
                assert d.dataset_status == expected_status

    def _test_change_status(self, method_name, required_status=None):
        with empty_analysis_driver_procs():
            d = RunDataset(
                name='test',
                path=os.path.join(self.base_dir, 'test'),
                use_int_dir=False
            )
            assert d.dataset_status == DATASET_NEW

        if not required_status:
            required_status = DATASET_READY
        with fake_analysis_driver_proc(d, required_status) as mocked_instance:
            method = d.__getattribute__(method_name)
            method()
            return d, mocked_instance

    def test_start(self):
        d, mocked_instance = self._test_change_status('start', DATASET_NEW)
        mocked_instance.assert_any_call('POST', cfg.query('rest_api', 'url') + '/analysis_driver_procs/', json={'status': DATASET_PROCESSING, 'end_date': d._now(), 'dataset_name': 'test', 'dataset_type': 'run'})
        mocked_instance.assert_any_call('POST', cfg.query('rest_api', 'url') + '/runs', json={'analysis_driver_procs': ['run_test_' + d._now()]})

    def test_abort(self):
        d, mocked_instance = self._test_change_status('abort')

    def test_fail(self):
        d, mocked_instance = self._test_change_status('fail', DATASET_PROCESSING)

    def test_reset(self):
        d, mocked_instance = self._test_change_status('reset')

    def test_succeed(self):
        d, mocked_instance = self._test_change_status('succeed', DATASET_PROCESSING)


class TestRunScanner(TestAnalysisDriver):
    @property
    def triggerignore(self):
        return os.path.join(self.scanner.input_dir, '.triggerignore')

    def setUp(self):
        self.base_dir = os.path.join(self.assets_path, 'dataset_scanner')
        seed_directories(self.base_dir)

        self.scanner = RunScanner({'lock_file_dir': self.base_dir, 'input_dir': self.base_dir})

        with open(self.triggerignore, 'w') as f:
            for d in ['ignored_dataset\n']:
                f.write(d)

    def tearDown(self):
        clean(self.base_dir)

    @patch('requests.request', new=fake_scan_request)
    def test_scan_datasets(self):
        datasets = self.scanner.scan_datasets()
        print('datasets: ', datasets)

        for status, expected in (
            (DATASET_NEW, {'dataset_not_ready'}),
            (DATASET_READY, {'dataset_ready'}),
            (DATASET_PROCESSING, {'processing_dataset'}),
            (DATASET_PROCESSED_FAIL, {'failed_dataset'})
        ):
            self.compare_lists(observed=self._dataset_names(datasets, status), expected=expected)

    @patch('requests.request', new=fake_scan_request)
    def test_triggerignore(self):
        with open(self.triggerignore, 'r') as f:
            assert f.readlines() == ['ignored_dataset\n']
            assert os.path.isdir(os.path.join(self.base_dir, 'ignored_dataset'))
            assert'ignored_dataset' not in self._flatten(self.scanner.scan_datasets())

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
