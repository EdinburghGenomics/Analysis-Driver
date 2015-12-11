__author__ = 'mwham'
import json
import os
from glob import glob
from collections import defaultdict
from analysis_driver.report_generation import ELEMENT_NB_Q30_R1, ELEMENT_NB_Q30_R2, ELEMENT_RUN_NAME
from analysis_driver.app_logging import get_logger


app_logger = get_logger('scanner')

DATASET_NEW = 'new'
DATASET_READY = 'ready'
DATASET_FORCE_READY = 'force_ready'
DATASET_PROCESSING = 'processing'
DATASET_PROCESSED_SUCCESS = 'finished'
DATASET_PROCESSED_FAIL = 'failed'
DATASET_ABORTED = 'aborted'

STATUS_VISIBLE = [DATASET_NEW, DATASET_READY, DATASET_FORCE_READY, DATASET_PROCESSING]
STATUS_HIDDEN = [DATASET_PROCESSED_SUCCESS, DATASET_PROCESSED_FAIL, DATASET_ABORTED]


class Dataset:
    def __init__(self, name, path, lock_file_dir):
        self.name = name
        self.path = path
        self.lock_file_dir = lock_file_dir
        self._stages = None
        self._pid = None

    @property
    def dataset_status(self):
        dataset_lock_files = glob(self._lock_file('*'))
        assert len(dataset_lock_files) < 2

        if dataset_lock_files:
            lf_status = dataset_lock_files[0].split('.')[-1]
        else:
            lf_status = DATASET_NEW

        if self._is_ready() and lf_status == DATASET_NEW:
            return DATASET_READY
        else:
            return lf_status

    @property
    def _is_ready(self):
        raise NotImplementedError

    def start(self):
        assert self.dataset_status in (DATASET_READY, DATASET_FORCE_READY, DATASET_NEW)
        self._change_status(DATASET_PROCESSING)
        self.set_pid()

    def succeed(self):
        assert self.dataset_status == DATASET_PROCESSING
        self.clear_pid()
        self._change_status(DATASET_PROCESSED_SUCCESS)

    def fail(self):
        assert self.dataset_status == DATASET_PROCESSING
        self._clear_stage()
        self.clear_pid()
        self._change_status(DATASET_PROCESSED_FAIL)

    def abort(self):
        self._clear_stage()
        self.clear_pid()
        self._change_status(DATASET_ABORTED)

    def reset(self):
        self._clear_stage()
        self.clear_pid()
        self._rm(*glob(self._lock_file('*')))

    def reset_status(self):
        self._rm(*glob(self._lock_file('*')))

    def _change_status(self, status):
        self.reset_status()
        self._touch(self._lock_file(status))

    def _lock_file(self, status):
        return os.path.join(
            self.lock_file_dir,
            '.' + self.name + '.' + status
        )

    def _stage_file(self, stage):
        return os.path.join(
            self.lock_file_dir,
            '.stage_' + self.name + '.' + stage
        )

    def _clear_stage(self):
        self._rm(*glob(self._stage_file('*')))

    def add_stage(self, stage):
        self._touch(self._stage_file(stage))

    def remove_stage(self, stage):
        self._rm(self._stage_file(stage))

    def _pid_file(self, pid):
        return os.path.join(self.lock_file_dir, '.pid_' + self.name + '.' + pid)

    def set_pid(self):
        self._touch(self._pid_file(str(os.getpid())))

    def clear_pid(self):
        self._rm(*glob(self._pid_file('*')))

    @property
    def stages(self):
        if self._stages is None:
            stage_files = glob(self._stage_file('*'))
            self._stages = [sf.split('.')[-1] for sf in stage_files]
        return self._stages

    @property
    def pid(self):
        if self._pid is None:
            pid_files = glob(self._pid_file('*'))
            if len(pid_files) == 1:
                self._pid = [sf.split('.')[-1] for sf in pid_files][0]
            else:
                self._pid = None
        return self._pid

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

    def _is_ready(self):
        return self.rta_complete()

    def rta_complete(self):
        return os.path.isfile(os.path.join(self.path, 'RTAComplete.txt'))


class SampleDataset(Dataset):
    type = 'sample'

    def __init__(self, name, path, lock_file_dir, data_threshold=None):
        super().__init__(name, path, lock_file_dir)
        self.data_threshold = data_threshold
        self.run_elements = self._read_data()

    def force(self):
        self._clear_stage()
        self.clear_pid()
        self._change_status(DATASET_FORCE_READY)

    def _read_data(self):
        with open(self.path, 'r') as open_file:
            return json.load(open_file)

    def _amount_data(self):
        return sum(
            [
                int(r.get(ELEMENT_NB_Q30_R1)) + int(r.get(ELEMENT_NB_Q30_R2))
                for r in self.run_elements.values()
            ]
        )

    def _runs(self):
        return set([r.get(ELEMENT_RUN_NAME) for r in self.run_elements.values()])

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
                        ignorables.extend(glob(os.path.join(self.input_dir, p.rstrip('\n'))))
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
    def _get_dataset(self, dataset_path):
        return RunDataset(
            name=os.path.basename(dataset_path),
            path=dataset_path,
            lock_file_dir=self.lock_file_dir
        )


class SampleScanner(DatasetScanner):
    def __init__(self, cfg):
        super().__init__(cfg)
        self.lock_file_dir = cfg.get('lock_file_dir', cfg['metadata_input_dir'])
        self.input_dir = cfg.get('metadata_input_dir')
        self.data_threshold = cfg.get('data_threshold')

    def _get_dataset(self, dataset_path):
        return SampleDataset(
            name=os.path.basename(dataset_path),
            path=dataset_path,
            lock_file_dir=self.lock_file_dir,
            data_threshold=self.data_threshold
        )
