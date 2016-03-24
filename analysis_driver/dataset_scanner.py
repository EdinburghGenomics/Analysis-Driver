import os
from datetime import datetime
from collections import defaultdict
from analysis_driver.exceptions import AnalysisDriverError
from analysis_driver import rest_communication
from analysis_driver.app_logging import get_logger
from analysis_driver.clarity import get_expected_yield_for_sample
from analysis_driver.constants import DATASET_NEW, DATASET_READY, DATASET_FORCE_READY, DATASET_PROCESSING,\
    DATASET_PROCESSED_SUCCESS, DATASET_PROCESSED_FAIL, DATASET_ABORTED, DATASET_REPROCESS, ELEMENT_RUN_NAME,\
    ELEMENT_NB_Q30_R2_CLEANED, ELEMENT_NB_Q30_R1_CLEANED

app_logger = get_logger('scanner')

STATUS_VISIBLE = [DATASET_NEW, DATASET_READY, DATASET_FORCE_READY, DATASET_PROCESSING]
STATUS_HIDDEN = [DATASET_PROCESSED_SUCCESS, DATASET_PROCESSED_FAIL, DATASET_ABORTED]


class Dataset:
    type = None
    endpoint = None
    id_field = None

    def __init__(self, name):
        self.name = name
        most_recent_proc = self._most_recent_proc()
        self.pid = most_recent_proc.get('pid')
        self.proc_id = most_recent_proc.get('proc_id', '_'.join((self.type, self.name)))

    def _most_recent_proc(self):
        procs = rest_communication.get_documents(
            'analysis_driver_procs',
            where={'dataset_type': self.type, 'dataset_name': self.name},
            sort='-_created'
        )
        if procs:
            return procs[0]
        else:
            return {}

    def _create_process(self, status, end_date=None):
        proc = {
            'proc_id': self.proc_id,
            'dataset_type': self.type,
            'dataset_name': self.name,
            'status': status
        }
        if end_date:
            proc['end_date'] = end_date
        if self.pid:
            proc['pid'] = self.pid
        rest_communication.post_entry('analysis_driver_procs', [proc])
        dataset = {self.id_field: self.name, 'analysis_driver_procs': [self.proc_id]}
        rest_communication.post_or_patch(
            self.endpoint,
            [dataset],
            id_field=self.id_field,
            update_lists=['analysis_driver_procs']
        )
        return proc

    @property
    def dataset_status(self):
        most_recent_proc = self._most_recent_proc()
        db_proc_status = most_recent_proc.get('status')
        if db_proc_status in (DATASET_REPROCESS, None):
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
        assert self.dataset_status in (DATASET_READY, DATASET_FORCE_READY, DATASET_NEW)
        self.pid = os.getpid()
        # sleep(1.1)
        start_time = self._now()
        self.proc_id = '_'.join((self.type, self.name, start_time))
        # proc_id is now different, so _change_status should create a new analysis_driver_proc and register it
        # to an appropriate endpoint
        self._change_status(DATASET_PROCESSING, finish=False)

    def succeed(self):
        assert self.dataset_status == DATASET_PROCESSING
        self._change_status(DATASET_PROCESSED_SUCCESS)

    def fail(self):
        assert self.dataset_status == DATASET_PROCESSING
        self._change_status(DATASET_PROCESSED_FAIL)

    def abort(self):
        self._change_status(DATASET_ABORTED)

    def reset(self):
        new_content = {'proc_id': self.proc_id, 'status': DATASET_REPROCESS}
        rest_communication.post_or_patch('analysis_driver_procs', [new_content], id_field='proc_id')

    def _change_status(self, status, finish=True):
        new_content = {
            'dataset_type': self.type,
            'dataset_name': self.name,
            'status': status
        }
        if finish:
            self.pid = None
            end_date = self._now()
            new_content['pid'] = self.pid
            new_content['end_date'] = end_date
        else:
            end_date = None

        patch_success = rest_communication.patch_entry(
            'analysis_driver_procs',
            new_content,
            'proc_id',
            self.proc_id
        )
        if not patch_success:
            self._create_process(status=status, end_date=end_date)

    def add_stage(self, stage_name):
        now = self._now()
        stages = self._most_recent_proc().get('stages', [])
        new_stage = {
            'date_started': now,
            'stage_name': stage_name
        }
        stages.append(new_stage)
        new_content = {'proc_id': self.proc_id, 'stages': stages}
        rest_communication.post_or_patch('analysis_driver_procs', [new_content], id_field='proc_id')

    def end_stage(self, stage_name, exit_status):
        stages = self._most_recent_proc().get('stages')
        for s in stages:
            if s['stage_name'] == stage_name:
                s['date_finished'] = self._now()
                s['exit_status'] = exit_status

        new_content = {'proc_id': self.proc_id, 'stages': stages}
        rest_communication.post_or_patch('analysis_driver_procs', [new_content], id_field='proc_id')

    @property
    def stages(self):
        proc = self._most_recent_proc()
        return [s['stage_name'] for s in proc.get('stages', []) if 'date_finished' not in s]

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
    endpoint = 'runs'
    id_field = 'run_id'

    def __init__(self, name, path, use_int_dir):
        super().__init__(name)
        self.path = path
        self.use_int_dir = use_int_dir

    def _is_ready(self):
        return self.rta_complete() or self.use_int_dir

    def rta_complete(self):
        return os.path.isfile(os.path.join(self.path, 'RTAComplete.txt'))


class SampleDataset(Dataset):
    type = 'sample'
    endpoint = 'samples'
    id_field = 'sample_id'

    def __init__(self, name):
        super().__init__(name)
        self.run_elements = self._read_data()
        self._data_threshold = None

    def force(self):
        self._change_status(DATASET_FORCE_READY, finish=False)

    def _read_data(self):
        return rest_communication.get_documents(
            'run_elements',
            where={'sample_id': self.name, 'useable': 'yes'}
        )

    def _amount_data(self):
        return sum(
            [
                int(r.get(ELEMENT_NB_Q30_R1_CLEANED, 0)) + int(r.get(ELEMENT_NB_Q30_R2_CLEANED, 0))
                for r in self.run_elements
            ]
        )

    def _runs(self):
        return sorted(set([r.get(ELEMENT_RUN_NAME) for r in self.run_elements]))

    @property
    def data_threshold(self):
        if self._data_threshold is None:
            self._data_threshold = get_expected_yield_for_sample(self.name)
        if not self._data_threshold:
            raise AnalysisDriverError('Could not find data threshold in LIMS for ' + self.name)
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
    def __init__(self, config):
        self.input_dir = config.get('input_dir')

    def scan_datasets(self):
        triggerignore = os.path.join(self.input_dir, '.triggerignore')

        ignorables = []
        if os.path.isfile(triggerignore):
            with open(triggerignore, 'r') as f:
                for p in f.readlines():
                    if not p.startswith('#'):
                        ignorables.append(p.rstrip('\n'))
        app_logger.debug('Ignoring %s datasets' % len(ignorables))

        n_datasets = 0
        datasets = defaultdict(list)
        for name in self._list_datasets():
            if name not in ignorables and not name.startswith('.'):
                d = self.get_dataset(name)
                datasets[d.dataset_status].append(d)
                n_datasets += 1
        app_logger.debug('Found %s datasets' % n_datasets)
        return datasets

    def _list_datasets(self):
        raise NotImplementedError

    def get_dataset(self, dataset_path):
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
    def __init__(self, config):
        super().__init__(config)
        self.use_int_dir = 'intermediate_dir' in config

    def _list_datasets(self):
        return os.listdir(self.input_dir)

    def get_dataset(self, name):
        dataset_path = os.path.join(self.input_dir, name)
        if os.path.exists(dataset_path):
            return RunDataset(
                name=os.path.basename(dataset_path),
                path=dataset_path,
                use_int_dir=self.use_int_dir
            )


class SampleScanner(DatasetScanner):
    def __init__(self, config):
        super().__init__(config)
        self.data_threshold = config.get('data_threshold')

    def _list_datasets(self, query=None):
        return [s['sample_id'] for s in rest_communication.get_documents('samples', depaginate=True)]

    def get_dataset(self, name):
        dataset_path = os.path.join(self.input_dir, name)
        return SampleDataset(name=os.path.basename(dataset_path))
