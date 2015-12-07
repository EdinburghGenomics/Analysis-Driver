import json
from analysis_driver.report_generation import ELEMENT_NB_Q30_R1, ELEMENT_NB_Q30_R2

__author__ = 'mwham'
import os
from glob import glob
from collections import defaultdict
from analysis_driver.app_logging import get_logger

app_logger = get_logger('scanner')

DATASET_NEW = 'new'
DATASET_READY = 'ready'
DATASET_PROCESSING = 'processing'
DATASET_PROCESSED_SUCCESS = 'finished'
DATASET_PROCESSED_FAIL = 'failed'
DATASET_ABORTED = 'aborted'

STATUS_VISIBLE=[DATASET_NEW, DATASET_READY, DATASET_PROCESSING]
STATUS_HIDEN=[DATASET_PROCESSED_SUCCESS, DATASET_PROCESSED_FAIL, DATASET_ABORTED]

class Dataset:
    def __init__(self, name, path, lock_file_dir, data_threshold=None):
        self.name = name
        self.path = path
        self.lock_file_dir = lock_file_dir
        self.data_threshold = data_threshold


    @property
    def dataset_status(self):
        raise NotImplementedError("Function not implemented in DatasetScanner")

    def start(self):
        assert self.dataset_status==DATASET_READY or self.dataset_status==DATASET_NEW
        self._change_status(DATASET_PROCESSING)
        self.set_pid()

    def succeed(self):
        assert self.dataset_status==DATASET_PROCESSING
        self.clear_pid()
        self._change_status(DATASET_PROCESSED_SUCCESS)

    def fail(self):
        assert self.dataset_status==DATASET_PROCESSING
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
        return os.path.join(
            self.lock_file_dir,
        '.pid_'+ self.name + '.' + pid)

    def set_pid(self):
        self._touch(self._pid_file(str(os.getpid())))

    def clear_pid(self):
        self._rm(*glob(self._pid_file('*')))


    @property
    def stages(self):
        if not hasattr(self,'_stages'):
            stage_files = glob(self._stage_file('*'))
            self._stages = [sf.split('.')[-1] for sf in stage_files]
        return self._stages

    @property
    def pid(self):
        if not hasattr(self,'_pid'):
            pid_files = glob(self._pid_file('*'))
            if len(pid_files) == 1 :
                self._pid = [sf.split('.')[-1] for sf in pid_files][0]
            else:
                self._pid=None
        return self._pid

    def _touch(self, file):
        open(file, 'w').close()

    def _rm(self, *files):
        for f in files:
            if os.path.isfile(f):
                os.remove(f)

    def __str__(self):
        out = [self.name]
        if self.pid:
            out.append('(%s)'%self.pid)
        if self.stages:
            out.append('-- %s'%(', '.join(self.stages)))
        return ' '.join(out)

    __repr__ = __str__

class RunDataset(Dataset):
    type = 'run'
    @property
    def dataset_status(self):
        dataset_lock_files = glob(self._lock_file('*'))
        assert len(dataset_lock_files) < 2

        if dataset_lock_files:
            lf_status = dataset_lock_files[0].split('.')[-1]
        else:
            lf_status = DATASET_NEW

        rta_complete = self._rta_complete()
        if rta_complete and lf_status == DATASET_NEW:
            return DATASET_READY
        else:
            return lf_status

    def _rta_complete(self):
        return os.path.isfile(os.path.join(self.path, 'RTAComplete.txt'))


class SampleDataset(Dataset):
    type = 'sample'
    @property
    def dataset_status(self):
        dataset_lock_files = glob(self._lock_file('*'))
        assert len(dataset_lock_files) < 2

        if dataset_lock_files:
            lf_status = dataset_lock_files[0].split('.')[-1]
        else:
            lf_status = DATASET_NEW

        if lf_status == DATASET_NEW and self._is_ready():
            return DATASET_READY
        else:
            return lf_status

    def _read_json(self):
        with open(self.path) as open_file:
            self.run_elements = json.load(open_file)

    def _amount_data(self):
        if not hasattr(self,'run_elements'):
            self._read_json()
        return sum([int(r.get(ELEMENT_NB_Q30_R1)) + int(r.get(ELEMENT_NB_Q30_R2)) for r in self.run_elements.values()])

    def _is_ready(self):
        if self.data_threshold and int(self._amount_data()) > int(self.data_threshold):
            return True
        else:
            return False

    def __str__(self):
        return '%s  (%s / %s)'%(super().__str__(), self._amount_data(), self.data_threshold)

class DatasetScanner():
    def __init__(self, cfg):
        self.lock_file_dir = cfg.get('lock_file_dir', cfg['input_dir'])
        self.input_dir = cfg.get('input_dir')
        self.data_threshold = cfg.get('data_threshold', None)

    def scan_datasets(self):
        raise NotImplementedError()

    def report(self, all_datasets=False):
        datasets = self.scan_datasets()
        out = []
        out.append('dataset location: ' + self.input_dir)
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
        return out

    def get(self, dataset_name, dataset_class):
        directory = glob(os.path.join(self.input_dir, dataset_name))
        if directory:
            directory = directory[0]
            d = dataset_class(name=os.path.basename(directory),
                        path=directory,
                        lock_file_dir=self.lock_file_dir)
            return d
        return None

class RunScanner(DatasetScanner):

    def __init__(self, cfg):
        super().__init__(cfg)


    def report(self, all_datasets=False):
        out = ['========= Run Scanner report =========']
        out.extend(super().report(all_datasets=all_datasets))
        print('\n'.join(out))

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
        for directory in glob(os.path.join(self.input_dir, '*')):
            d = RunDataset(
                name=os.path.basename(directory),
                path=directory,
                lock_file_dir=self.lock_file_dir
            )
            if os.path.isdir(directory) and directory not in ignorables:
                datasets[d.dataset_status].append(d)
                n_datasets += 1
        app_logger.debug('Found %s datasets' % n_datasets)
        return datasets

    def get(self, dataset_name):
        return super().get(dataset_name, RunDataset)

class SampleScanner(DatasetScanner):

    def __init__(self, cfg):
        self.lock_file_dir = cfg.get('lock_file_dir', cfg['metadata_input_dir'])
        self.input_dir = cfg.get('metadata_input_dir')
        self.data_threshold = cfg.get('data_threshold', None)


    def report(self, all_datasets=False):
        out = ['========= Sample Scanner report =========']
        out.extend(super().report(all_datasets=all_datasets))
        print('\n'.join(out))

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
        for directory in glob(os.path.join(self.input_dir, '*')):
            if os.path.isfile(directory) and directory not in ignorables:
                d = SampleDataset(
                    name=os.path.basename(directory),
                    path=directory,
                    lock_file_dir=self.lock_file_dir,
                    data_threshold=self.data_threshold
                )
                datasets[d.dataset_status].append(d)
                n_datasets += 1
        app_logger.debug('Found %s datasets' % n_datasets)
        return datasets

    def get(self, dataset_name):
        return super().get(dataset_name, SampleDataset)
