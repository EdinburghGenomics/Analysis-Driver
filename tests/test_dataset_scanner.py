__author__ = 'tcezard'
import shutil
import pytest
from unittest.mock import patch, Mock, PropertyMock
import os
from tests.test_analysisdriver import TestAnalysisDriver
from analysis_driver import util
from analysis_driver.exceptions import AnalysisDriverError
from analysis_driver.dataset_scanner import DatasetScanner, RunScanner, SampleScanner, DATASET_NEW, DATASET_READY, DATASET_PROCESSING, \
    DATASET_PROCESSED_FAIL, DATASET_PROCESSED_SUCCESS, DATASET_ABORTED, DATASET_REPROCESS, DATASET_FORCE_READY,\
    Dataset, RunDataset, SampleDataset
from analysis_driver.config import default as cfg


directories_to_create = ('dataset_ready', 'dataset_not_ready', 'ignored_dataset')
ready_datasets = ('dataset_ready',)


def seed_directories(base_dir):
    for d in directories_to_create:
        os.makedirs(os.path.join(base_dir, d), exist_ok=True)
        if d in ready_datasets:
            _touch(os.path.join(base_dir, d, 'RTAComplete.txt'))


def clean(base_dir):
    shutil.rmtree(os.path.join(base_dir))


def _touch(path):
    open(path, 'w').close()


def api(endpoint):
    return cfg.query('rest_api', 'url') + '/' + endpoint


fake_analysis_driver_proc = {
    'dataset_type': 'a_type', 'dataset_name': 'a_name', 'proc_id': 'a_type_a_name', 'status': 'a_status'
}
fake_analysis_driver_proc_no_status = {
    'dataset_type': 'a_type', 'dataset_name': 'a_name', 'proc_id': 'a_type_a_name'
}
fake_analysis_driver_proc_stages = {
    'dataset_type': 'a_type',
    'dataset_name': 'a_name',
    'proc_id': 'a_type_a_name',
    'stages': [{'date_started': 'a_start_date', 'stage_name': 'a_stage'}]
}
fake_sample = {
    'library_id': 'a_library_id',
    'project_id': 'a_project_id',
    'sample_id': 'a_sample_id',
    'user_sample_id': 'a_user_sample_id',
}


class FakeRestResponse(Mock):
    def json(self):
        return self.content


patched_request = patch(
    'requests.request',
    return_value=FakeRestResponse(content={'_links': {}, 'data': [fake_analysis_driver_proc]})
)
patched_request_no_data = patch(
    'requests.request',
    return_value=FakeRestResponse(content={'_links': {}, 'data': []})
)

patched_post_or_patch = patch('analysis_driver.rest_communication.post_or_patch')
patched_change_status = patch('analysis_driver.dataset_scanner.Dataset._change_status')
patched_stages = patch('analysis_driver.dataset_scanner.Dataset.stages', new_callable=PropertyMock(return_value=['this', 'that', 'other']))
patched_expected_yield = patch(
    'analysis_driver.dataset_scanner.get_expected_yield_for_sample', return_value=1000000000
)
patched_get_fake_sample = patch(
    'analysis_driver.dataset_scanner.requests.get',
    return_value=FakeRestResponse(content={'_links': {}, 'data': [fake_sample]})
)


def patched_dataset_status(status=DATASET_NEW):
    return patch(
        'analysis_driver.dataset_scanner.Dataset.dataset_status',
        new_callable=PropertyMock(return_value=status)
    )


def patched_most_recent_proc(proc=fake_analysis_driver_proc):
    return patch(
        'analysis_driver.dataset_scanner.Dataset._most_recent_proc',
        return_value=proc
    )


def patched_get_docs(content):
    return patch('analysis_driver.rest_communication.get_documents', return_value=content)


class TestDataset(TestAnalysisDriver):
    def setUp(self):
        self.base_dir = os.path.join(self.assets_path, 'dataset_scanner')
        seed_directories(self.base_dir)
        self.setup_dataset()

    def tearDown(self):
        clean(self.base_dir)

    def test_most_recent_proc(self):
        expected_api_call = util.str_join(
            api('analysis_driver_procs'),
            '?where={"dataset_type":"',
            self.dataset.type,
            '","dataset_name":"test_dataset"}&sort=-_created'
        )
        with patched_request as mocked_instance:
            proc = self.dataset._most_recent_proc()
            assert proc == fake_analysis_driver_proc
            mocked_instance.assert_called_with('GET', expected_api_call)

        with patched_request_no_data as mocked_instance:
            proc = self.dataset._most_recent_proc()
            assert proc == {}
            mocked_instance.assert_called_with('GET', expected_api_call)

    def test_create_process(self):
        with patch('analysis_driver.rest_communication.post_entry', return_value=True) as mocked_patch:
            with patched_post_or_patch as mocked_post_or_patch:
                proc = self.dataset._create_process('a_status', end_date='an_end_date')

                assert proc == {
                    'proc_id': self.dataset.proc_id,
                    'dataset_type': self.dataset.type,
                    'dataset_name': self.dataset.name,
                    'status': 'a_status',
                    'end_date': 'an_end_date'
                }
                mocked_patch.assert_called_with('analysis_driver_procs', [proc])
                mocked_post_or_patch.assert_called_with(
                    self.dataset.endpoint,
                    [{self.dataset.id_field: 'test_dataset', 'analysis_driver_procs': ['a_type_a_name']}],
                    elem_key=self.dataset.id_field,
                    update_lists=['analysis_driver_procs']
                )

    def test_dataset_status(self):
        with patched_most_recent_proc() as mocked_instance:
            assert self.dataset.dataset_status == 'a_status'
            assert mocked_instance.call_count == 1
            assert mocked_instance.return_value == {
                'status': 'a_status',
                'proc_id': 'a_type_a_name',
                'dataset_type': 'a_type',
                'dataset_name': 'a_name'
            }

    @patch('analysis_driver.dataset_scanner.Dataset._create_process')
    def test_change_status(self, mocked_instance):
        with patch('analysis_driver.rest_communication.patch_entry', return_value=False) as mocked_patch:
            self.dataset._change_status('a_status', finish=True)
            mocked_patch.assert_called_with(
                'analysis_driver_procs',
                {
                    'status': 'a_status',
                    'dataset_type': self.dataset.type,
                    'dataset_name': 'test_dataset',
                    'end_date': self.dataset._now()
                },
                proc_id='a_type_a_name'
            )
            mocked_instance.assert_called_with(status='a_status', end_date=self.dataset._now())

    @patched_change_status
    def test_start(self, mocked_change_status):
        with patched_dataset_status():
            assert self.dataset.proc_id == 'a_type_a_name'
            self.dataset.start()
            assert self.dataset.proc_id == self.dataset.type + '_test_dataset_' + self.dataset._now()
            mocked_change_status.assert_called_with(DATASET_PROCESSING, finish=False)

        with patched_dataset_status(DATASET_ABORTED), pytest.raises(AssertionError):
            self.dataset.start()

    @patched_change_status
    def test_succeed(self, mocked_change_status):
        with patched_dataset_status(DATASET_PROCESSING):
            self.dataset.succeed()
            mocked_change_status.assert_called_with(DATASET_PROCESSED_SUCCESS)

        with patched_dataset_status(DATASET_NEW), pytest.raises(AssertionError):
            self.dataset.succeed()

    @patched_change_status
    def test_fail(self, mocked_change_status):
        with patched_dataset_status(DATASET_PROCESSING):
            self.dataset.fail()
            mocked_change_status.assert_called_with(DATASET_PROCESSED_FAIL)

        with patched_dataset_status(DATASET_NEW), pytest.raises(AssertionError):
            self.dataset.fail()

    @patched_change_status
    def test_abort(self, mocked_change_status):
        with patched_dataset_status():
            self.dataset.abort()
            mocked_change_status.assert_called_with(DATASET_ABORTED)

    @patched_post_or_patch
    def test_reset(self, mocked_instance):
        self.dataset.reset()
        mocked_instance.assert_called_with(
            'analysis_driver_procs',
            [{'proc_id': 'a_type_a_name', 'status': DATASET_REPROCESS}],
            elem_key='proc_id'
        )

    @patched_most_recent_proc()
    @patched_post_or_patch
    def test_add_stage(self, mocked_post_or_patch, mocked_most_recent_proc):
        now = self.dataset._now()
        self.dataset.add_stage('a_stage')
        mocked_post_or_patch.assert_called_with(
            'analysis_driver_procs',
            [
                {
                    'proc_id': 'a_type_a_name',
                    'stages': [{'date_started': now, 'stage_name': 'a_stage'}]
                }
            ],
            elem_key='proc_id'
        )

    @patched_most_recent_proc(fake_analysis_driver_proc_stages)
    @patched_post_or_patch
    def test_end_stage(self, mocked_post_or_patch, mocked_most_recent_proc):
        self.dataset.end_stage('a_stage', 0)
        mocked_post_or_patch.assert_called_with(
            'analysis_driver_procs',
            [
                {
                    'proc_id': 'a_type_a_name',
                    'stages': [
                        {
                            'date_started': 'a_start_date',
                            'stage_name': 'a_stage',
                            'exit_status': 0,
                            'date_finished': self.dataset._now()
                        }
                    ]
                }
            ],
            elem_key='proc_id'
        )

    def test_stages(self):
        fake_proc = {
            'dataset_type': 'a_type',
            'dataset_name': 'a_name',
            'proc_id': 'a_type_a_name',
            'status': 'a_status',
            'stages': [
                {
                    'stage_name': 'a_stage',
                    'date_started': '01_01_2016_12:00:00',
                    'date_finished': '01_01_2016_13:00:00'
                },
                {
                    'stage_name': 'another_stage',
                    'date_started': '03_01_2016_12:00:00'
                },
                {
                    'stage_name': 'yet_another_stage',
                    'date_started': '03_01_2016_13:00:00'
                }
            ]
        }
        with patched_most_recent_proc(fake_proc):
            assert self.dataset.stages == ['another_stage', 'yet_another_stage']

    @patched_stages
    def test_str(self, mocked_instance):
        assert str(self.dataset) == 'test_dataset -- this, that, other'

    def setup_dataset(self):
        with patched_request:
            self.dataset = Dataset('test_dataset', os.path.join(self.base_dir, 'test_dataset'))


class TestRunDataset(TestDataset):
    @patched_request
    def test_is_ready(self, mocked_instance):
        datasets = (
            ('dataset_ready', True),
            ('dataset_not_ready', False)
        )
        for d_name, rta_complete in datasets:
            d = RunDataset(d_name, os.path.join(self.base_dir, d_name), use_int_dir=False)
            assert d._is_ready() == rta_complete

    def test_dataset_status(self):
        super().test_dataset_status()
        with patched_most_recent_proc(fake_analysis_driver_proc_no_status):
            assert not self.dataset.rta_complete()
            assert self.dataset.dataset_status == DATASET_NEW
            os.mkdir(os.path.join(self.base_dir, self.dataset.name))
            _touch(os.path.join(self.base_dir, self.dataset.name, 'RTAComplete.txt'))
            assert self.dataset.rta_complete()
            assert self.dataset.dataset_status == DATASET_READY

    def setup_dataset(self):
        with patched_request:
            self.dataset = RunDataset(
                'test_dataset',
                os.path.join(self.base_dir, 'test_dataset'),
                use_int_dir=False
            )


class TestSampleDataset(TestDataset):
    @patched_change_status
    def test_force(self, mocked_instance):
        self.dataset.force()
        mocked_instance.assert_called_with(DATASET_FORCE_READY, finish=False)

    def test_amount_data(self):
        assert self.dataset._amount_data() == 480

    def test_runs(self):
        with patched_get_docs(self.dataset.run_elements):
            assert self.dataset._runs() == ['a_run_id', 'another_run_id']

    @patched_expected_yield
    def test_data_threshold(self, mocked_exp_yield):
        assert self.dataset.data_threshold == 1000000000
        mocked_exp_yield.assert_called_with('test_dataset')

    @patch('analysis_driver.dataset_scanner.get_expected_yield_for_sample', return_value=None)
    def test_no_data_threshold(self, mocked_exp_yield):
        with pytest.raises(AnalysisDriverError) as e:
            self.dataset.data_threshold
        assert 'Could not find data threshold in LIMS' in str(e)
        mocked_exp_yield.assert_called_with('test_dataset')

    @patched_expected_yield
    def test_is_ready(self, mocked_instance):
        assert not self.dataset._is_ready()
        self.dataset.run_elements = [
            {
                'clean_bases_r1': 1300000000,
                'clean_q30_bases_r1': 1200000000,
                'clean_bases_r2': 1250000000,
                'clean_q30_bases_r2': 1150000000
            },
            {
                'clean_bases_r1': 1200000000,
                'clean_q30_bases_r1': 1100000000,
                'clean_bases_r2': 1450000000,
                'clean_q30_bases_r2': 1350000000
            }
        ]
        assert self.dataset._is_ready()
        assert mocked_instance.call_count == 1  # even after 2 calls to data_threshold

    @patched_stages
    def test_str(self, mocked_instance):
        expected_str = 'test_dataset -- this, that, other  (480 / 1000000000  from a_run_id, another_run_id) '
        with patched_get_docs(self.dataset.run_elements), patched_expected_yield:
            print(expected_str)
            print(str(self.dataset))
            assert str(self.dataset) == expected_str

    def setup_dataset(self):
        with patched_request:
            self.dataset = SampleDataset('test_dataset', os.path.join(self.base_dir, 'test_dataset'))
            self.dataset.run_elements = [
                {
                    'run_element_id': 'a_run_element_id',
                    'run_id': 'a_run_id',
                    'lane': 1,
                    'barcode': 'TACGTACG',
                    'project_id': 'a_project_id',
                    'library_id': 'a_library_id',
                    'sample_id': 'a_sample_id',
                    'total_reads': 1337,
                    'passing_filter_reads': 1000,
                    'pc_reads_in_lane': 80,
                    'clean_bases_r1': 130,
                    'clean_q30_bases_r1': 120,
                    'clean_bases_r2': 125,
                    'clean_q30_bases_r2': 115
                },
                {
                    'run_element_id': 'a_run_element_id',
                    'run_id': 'another_run_id',
                    'lane': 1,
                    'barcode': 'ATGCATGC',
                    'project_id': 'a_project_id',
                    'library_id': 'a_library_id',
                    'sample_id': 'a_sample_id',
                    'total_reads': 1338,
                    'passing_filter_reads': 980,
                    'pc_reads_in_lane': 79,
                    'clean_bases_r1': 120,
                    'clean_q30_bases_r1': 110,
                    'clean_bases_r2': 145,
                    'clean_q30_bases_r2': 135
                }
            ]


class TestScanner(TestAnalysisDriver):
    @property
    def triggerignore(self):
        return os.path.join(self.scanner.input_dir, '.triggerignore')

    def setUp(self):
        self.base_dir = os.path.join(self.assets_path, 'dataset_scanner')
        seed_directories(self.base_dir)
        self._setup_scanner()
        assert self.scanner.input_dir == self.base_dir

        with open(self.triggerignore, 'w') as f:
            for d in ['ignored_dataset\n']:
                f.write(d)

    def tearDown(self):
        clean(self.base_dir)

    @patched_most_recent_proc({})
    def test_scan_datasets(self, mocked_instance):
        with pytest.raises(NotImplementedError):
            self.scanner.scan_datasets()

    def test_report(self):
        dsets = {'new': ['this'], 'ready': ['that', 'other'], 'failed': ['another']}
        captured_stdout = []
        with patch('analysis_driver.dataset_scanner.DatasetScanner.scan_datasets', return_value=dsets):
            with patch('builtins.print', new=captured_stdout.append):
                self.scanner.report(all_datasets=True)

        expected = [
            '========= ' + self.scanner.__class__.__name__ + ' report =========',
            'dataset location: ' + self.base_dir,
            '=== new ===',
            'this',
            '=== ready ===',
            'that',
            'other',
            '=== failed ===',
            'another',
            '__________________________________________'
        ]

        observed = captured_stdout[0].split('\n')
        print(observed)
        print(expected)
        assert observed == expected

    def _fake_get_dataset(self, name):
        with patched_most_recent_proc():
            return Dataset(os.path.basename(name), os.path.join(self.scanner.input_dir, name))

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

    def _setup_scanner(self):
        self.scanner = DatasetScanner({'input_dir': self.base_dir})


class TestRunScanner(TestScanner):
    @patched_request
    def test_get_dataset(self, mocked_instance):
        test_dataset_path = os.path.join(self.base_dir, 'test_dataset')
        os.mkdir(test_dataset_path)
        observed = self.scanner.get_dataset(test_dataset_path)
        expected = RunDataset('test_dataset', test_dataset_path, self.scanner.use_int_dir)
        for a in ('name', 'path', 'use_int_dir'):
            assert observed.__getattribute__(a) == expected.__getattribute__(a)

    def test_list_datasets(self):
        self.compare_lists(self.scanner._list_datasets(), os.listdir(self.base_dir))

    @patched_most_recent_proc({})
    def test_scan_datasets(self, mocked_instance):
        datasets = self.scanner.scan_datasets()
        print('datasets: ', datasets)

        for status, expected in (
            (DATASET_NEW, {'dataset_not_ready'}),
            (DATASET_READY, {'dataset_ready'}),
        ):
            self.compare_lists(observed=self._dataset_names(datasets, status), expected=expected)

    @patched_most_recent_proc()
    def test_triggerignore(self, mocked_instance):
        with open(self.triggerignore, 'r') as f:
            assert f.readlines() == ['ignored_dataset\n']
            assert os.path.isdir(os.path.join(self.base_dir, 'ignored_dataset'))
            assert 'ignored_dataset' not in self._flatten(self.scanner.scan_datasets())

    def _setup_scanner(self):
        self.scanner = RunScanner({'input_dir': self.base_dir})


class TestSampleScanner(TestScanner):
    def _setup_scanner(self):
        self.scanner = SampleScanner({'input_dir': self.base_dir})

    @patched_request
    def test_get_dataset(self, mocked_instance):
        with patched_expected_yield:
            observed = self.scanner.get_dataset(os.path.join(self.base_dir, 'test_dataset'))
            expected = SampleDataset('test_dataset', os.path.join(self.base_dir, 'test_dataset'))
            for a in ('name', 'path', 'data_threshold'):
                assert observed.__getattribute__(a) == expected.__getattribute__(a)

    @patched_get_fake_sample
    def test_list_datasets(self, mocked_instance):
        assert self.scanner._list_datasets() == ['a_sample_id']
        mocked_instance.assert_called_with(api('samples'))

    @patched_request
    def test_scan_datasets(self, mocked_instance):
        with patched_most_recent_proc(fake_analysis_driver_proc_no_status):
            statuses = (
                DATASET_ABORTED, DATASET_NEW, DATASET_READY, DATASET_FORCE_READY, DATASET_PROCESSED_FAIL,
                DATASET_PROCESSED_SUCCESS, DATASET_PROCESSING, DATASET_REPROCESS
            )
            for status in statuses:
                print(status)
                with patched_dataset_status(status), patched_get_fake_sample:
                    datasets = self.scanner.scan_datasets()
                    assert list(datasets.keys()) == [status]
                    self.compare_lists([d.name for d in datasets[status]], ['a_sample_id'])
