__author__ = 'mwham'
import os
import requests
from datetime import datetime
from time import sleep
from collections import defaultdict
from analysis_driver.config import default as cfg
from analysis_driver.report_generation import rest_communication, ELEMENT_NB_Q30_R1, ELEMENT_NB_Q30_R2,\
    ELEMENT_RUN_NAME
from analysis_driver.app_logging import get_logger
from analysis_driver.clarity import get_expected_yield_for_sample


app_logger = get_logger('scanner')

DATASET_NEW = 'new'
DATASET_READY = 'ready'
DATASET_FORCE_READY = 'force_ready'
DATASET_PROCESSING = 'processing'
DATASET_PROCESSED_SUCCESS = 'finished'
DATASET_PROCESSED_FAIL = 'failed'
DATASET_ABORTED = 'aborted'
DATASET_REPROCESS = 'reprocess'

STATUS_VISIBLE = [DATASET_NEW, DATASET_READY, DATASET_FORCE_READY, DATASET_PROCESSING]
STATUS_HIDDEN = [DATASET_PROCESSED_SUCCESS, DATASET_PROCESSED_FAIL, DATASET_ABORTED]


class Dataset:
    type = None

    def __init__(self, name, path, lock_file_dir):
        self.name = name
        self.path = path
        self.lock_file_dir = lock_file_dir
        self._stages = None
        self.pid = None
        self.proc_id = self._most_recent_proc().get('proc_id', '_'.join((self.type, self.name)))

    def _most_recent_proc(self):
        # TODO: add embedding, sort, etc. support into rest_communication
        query_url = ''.join(
            (
                cfg.query('rest_api', 'url').rstrip('/'),
                '/analysis_driver_procs?where={"dataset_type":"',
                self.type,
                '","dataset_name":"',
                self.name,
                '"}&sort=-_created'
            )
        )
        procs = requests.request('GET', query_url).json()['data']
        if procs:
            return procs[0]
        else:
            return {}

    def _create_process(self, **extra_params):
        proc = {
            'proc_id': self.proc_id,
            'dataset_type': self.type,
            'dataset_name': self.name,
        }
        proc.update(extra_params)
        rest_communication.post_or_patch('analysis_driver_procs', [proc], elem_key='proc_id')
        name_key = self.type + '_id'
        dataset = {name_key: self.name, 'analysis_driver_procs': [self.proc_id]}
        rest_communication.post_or_patch(
            self.type + 's',
            [dataset],
            elem_key=name_key,
            update_lists=['analysis_driver_procs']
        )
        return proc

    @property
    def dataset_status(self):
        most_recent_proc = self._most_recent_proc()
        db_proc_status = most_recent_proc.get('status')
        if not db_proc_status:
            if self._is_ready():
                return DATASET_READY
            else:
                return DATASET_NEW
        else:
            return db_proc_status

    @property
    def _is_ready(self):
        raise NotImplementedError

    @staticmethod
    def _now():
        return datetime.utcnow().strftime('%d_%m_%Y_%H:%M:%S')

    def start(self):
        assert self.dataset_status in (DATASET_READY, DATASET_FORCE_READY, DATASET_NEW, DATASET_REPROCESS)
        self.pid = os.getpid()
        sleep(1.1)
        start_time = self._now()
        self.proc_id = '_'.join((self.type, self.name, start_time))
        # proc_id is now different, so _change_status should create a new analysis_driver_proc and register it
        # to an appropriate endpoint
        self._change_status(DATASET_PROCESSING)

    def succeed(self):
        self._change_status(DATASET_PROCESSED_SUCCESS, finish=True)

    def fail(self):
        assert self.dataset_status == DATASET_PROCESSING
        self._change_status(DATASET_PROCESSED_FAIL, finish=True)

    def abort(self):
        self._change_status(DATASET_ABORTED, finish=True)

    def reset(self):
        new_content = {'proc_id': self.proc_id, 'status': DATASET_REPROCESS}
        rest_communication.post_or_patch('analysis_driver_procs', [new_content], elem_key='proc_id')

    def _change_status(self, status, finish=False):
        now = self._now()
        new_content = {
            'proc_id': self.proc_id,
            'dataset_type': self.type,
            'dataset_name': self.name,
            'status': status
        }
        if finish:
            new_content['end_date'] = now

        patch_success = rest_communication.patch_entry(
            cfg.query('rest_api', 'url').rstrip('/') + '/analysis_driver_procs',
            new_content,
            proc_id=self.proc_id
        )
        if not patch_success:
            self._create_process(end_date=now, status=status)

    def add_stage(self, stage):
        new_content = {'proc_id': self.proc_id, 'stages': [stage]}
        rest_communication.post_or_patch('analysis_driver_procs', [new_content], update_lists=['stages'])

    def remove_stage(self, stage):
        pass

    @property
    def stages(self):
        if self._stages is None:
            proc = self._most_recent_proc()
            self._stages = [s['stage_name'] for s in proc.get('stages', []) if 'date_finished' in s]
        return self._stages

    @staticmethod
    def _touch(file):
        open(file, 'w').close()

    @staticmethod
    def _rm(*files):
        for f in files:
            if os.path.isfile(f):
                os.remove(f)

    def __str__(self):
        out = [self.name]
        if self.pid:
            out.append('(%s)' % self.pid)
        if self.stages:
            out.append('-- %s' % (', '.join(self.stages)))
        return ' '.join(out)

    __repr__ = __str__


class RunDataset(Dataset):
    type = 'run'

    def __init__(self, name, path, lock_file_dir, use_int_dir):
        super().__init__(name, path, lock_file_dir)
        self.use_int_dir = use_int_dir

    def _is_ready(self):
        return self.rta_complete() or self.use_int_dir

    def rta_complete(self):
        return os.path.isfile(os.path.join(self.path, 'RTAComplete.txt'))


class SampleDataset(Dataset):
    type = 'sample'

    def __init__(self, name, path, lock_file_dir, data_threshold=None):
        super().__init__(name, path, lock_file_dir)
        self.default_data_threshold = data_threshold
        self.run_elements = self._read_data()

    def force(self):
        self._change_status(DATASET_FORCE_READY)

    def _read_data(self):
        return rest_communication.get_documents(
            cfg['rest_api']['url'].rstrip('/') + '/run_elements',
            sample_id=self.name
        )

    def _amount_data(self):
        return sum(
            [
                int(r.get(ELEMENT_NB_Q30_R1)) + int(r.get(ELEMENT_NB_Q30_R2))
                for r in self.run_elements
            ]
        )

    def _runs(self):
        return set([r.get(ELEMENT_RUN_NAME) for r in self.run_elements])

    @property
    def data_threshold(self):
        if not hasattr(self, '_data_threshold'):
            self._data_threshold = get_expected_yield_for_sample(self.name)
        if not self._data_threshold:
            self._data_threshold = self.default_data_threshold
        return self._data_threshold

    def _is_ready(self):
        return self.data_threshold and int(self._amount_data()) > int(self.data_threshold)

    def __str__(self):
        return '%s  (%s / %s  from %s) ' % (
            super().__str__(),
            self._amount_data(),
            self.data_threshold,
            ', '.join(self._runs())
        )


class DatasetScanner:
    def __init__(self, cfg):
        self.lock_file_dir = cfg.get('lock_file_dir', cfg['input_dir'])
        self.input_dir = cfg.get('input_dir')

    def scan_datasets(self):
        triggerignore = os.path.join(self.lock_file_dir, '.triggerignore')

        ignorables = []
        if os.path.isfile(triggerignore):
            with open(triggerignore, 'r') as f:
                for p in f.readlines():
                    if not p.startswith('#'):
                        ignorables.append(p.rstrip('\n'))
        app_logger.debug('Ignoring %s datasets' % len(ignorables))

        n_datasets = 0
        datasets = defaultdict(list)
        for name in os.listdir(self.input_dir):
            if name not in ignorables and not name.startswith('.'):
                d = self.get(os.path.join(self.input_dir, name))
                datasets[d.dataset_status].append(d)
                n_datasets += 1
        app_logger.debug('Found %s datasets' % n_datasets)
        return datasets

    def get(self, name):
        dataset_path = os.path.join(self.input_dir, name)
        if os.path.exists(dataset_path):
            return self._get_dataset(dataset_path)

    def _get_dataset(self, dataset_path):
        raise NotImplementedError

    def report(self, all_datasets=False):
        out = [
            '========= %s report =========' % self.__class__.__name__,
            'dataset location: ' + self.input_dir
        ]
        datasets = self.scan_datasets()
        for status in STATUS_VISIBLE:
            ds = datasets.pop(status, [])
            if ds:
                out.append('=== ' + status + ' ===')
                out.append('\n'.join((str(d) for d in ds)))

        if any((datasets[s] for s in datasets)):
            if all_datasets:
                for status in sorted(datasets):
                    out.append('=== ' + status + ' ===')
                    out.append('\n'.join((str(d) for d in datasets[status])))
            else:
                out.append('=== other datasets ===')
                out.append('\n'.join(('other datasets present', 'use --report-all to show')))

        out.append('_' * 42)
        print('\n'.join(out))


class RunScanner(DatasetScanner):
    def __init__(self, cfg):
        super().__init__(cfg)
        self.use_int_dir = 'intermediate_dir' in cfg

    def _get_dataset(self, dataset_path):
        return RunDataset(
            name=os.path.basename(dataset_path),
            path=dataset_path,
            lock_file_dir=self.lock_file_dir,
            use_int_dir=self.use_int_dir
        )


class SampleScanner(DatasetScanner):
    def __init__(self, cfg):
        super().__init__(cfg)
        self.lock_file_dir = cfg.get('lock_file_dir', cfg['metadata_input_dir'])
        self.input_dir = cfg.get('metadata_input_dir')  # override input_dir
        self.data_threshold = cfg.get('data_threshold')

    def _get_dataset(self, dataset_path):
        return SampleDataset(
            name=os.path.basename(dataset_path),
            path=dataset_path,
            lock_file_dir=self.lock_file_dir,
            data_threshold=self.data_threshold
        )
