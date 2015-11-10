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
    def __init__(self, name, path, lock_file_dir):
        self.name = name
        self.path = path
        self.lock_file_dir = lock_file_dir

    @property
    def dataset_status(self):
        raise NotImplementedError("Function not implemented in DatasetScanner")

    def start(self):
        assert self.dataset_status==DATASET_READY
        self._change_status(DATASET_PROCESSING)

    def succeed(self):
        assert self.dataset_status==DATASET_PROCESSING
        self._change_status(DATASET_PROCESSED_SUCCESS)

    def fail(self):
        assert self.dataset_status==DATASET_PROCESSING
        self._change_status(DATASET_PROCESSED_FAIL)

    def abort(self):
        self._change_status( DATASET_ABORTED)

    def reset(self):
        self._rm(*glob(self._lock_file('*')))

    def _change_status(self, status):
        self.reset()
        self._touch(self._lock_file(status))

    def _lock_file(self, status):
        return os.path.join(
            self.lock_file_dir,
            '.' + self.name + '.' + status
        )

    def _touch(self, file):
        open(file, 'w').close()

    def _rm(self, *files):
        for f in files:
            if os.path.isfile(f):
                os.remove(f)

    def __str__(self):
        return self.name

    __repr__ = __str__

class RunDataset(Dataset):
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



class DatasetScanner():
    def __init__(self, cfg):
        self.lock_file_dir = cfg.get('lock_file_dir', cfg['input_dir'])
        self.input_dir = cfg.get('input_dir')

    def scan_datasets(self, dataset_class=Dataset):
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
            d = dataset_class(name=os.path.basename(directory),
                        path=directory,
                        lock_file_dir=self.lock_file_dir)
            if os.path.isdir(directory) and directory not in ignorables:
                datasets[d.dataset_status].append(d)
                n_datasets += 1
        app_logger.debug('Found %s datasets' % n_datasets)
        return datasets

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
        datasets = self.scan_datasets()
        out = []
        out.append('========= Run Scanner report =========')
        out.append('dataset location: ' + self.input_dir)
        for status in STATUS_VISIBLE:
            ds = datasets.pop(status, [])
            if ds:
                out.append('=== ' + status + ' ===')
                out.append('\n'.join((os.path.basename(d) for d in ds)))

        if any((datasets[s] for s in datasets)):
            if all_datasets:
                for status in sorted(datasets):
                    out.append('=== ' + status + ' ===')
                    out.append('\n'.join((os.path.basename(x) for x in datasets[status])))
            else:
                out.append('=== other datasets ===')
                out.append('\n'.join(('other datasets present', 'use --report-all to show')))

        out.append('_' * 42)
        print('\n'.join(out))

    def scan_datasets(self):
        return super().scan_datasets(RunDataset)

    def get(self, dataset_name):
        return super().get(dataset_name, RunDataset)

class SampleScanner(DatasetScanner):

    def __init__(self, cfg):
        super().__init__(cfg)


    def report(self, all_datasets=False):
        datasets = self.scan_datasets()
        out = []
        out.append('========= Sample Scanner report =========')
        out.append('dataset location: ' + self.input_dir)
        for status in ('new', 'new, rta complete', 'transferring', 'transferring, rta complete', 'active'):
            ds = datasets.pop(status, [])
            if ds:
                out.append('=== ' + status + ' ===')
                out.append('\n'.join((os.path.basename(d) for d in ds)))

        if any((datasets[s] for s in datasets)):
            if all_datasets:
                for status in sorted(datasets):
                    out.append('=== ' + status + ' ===')
                    out.append('\n'.join((os.path.basename(x) for x in datasets[status])))
            else:
                out.append('=== other datasets ===')
                out.append('\n'.join(('other datasets present', 'use --report-all to show')))

        out.append('_' * 42)

    def scan_datasets(self):
        return super().scan_datasets(SampleDataset)