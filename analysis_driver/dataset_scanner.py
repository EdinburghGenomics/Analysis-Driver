import os
from datetime import datetime
from collections import defaultdict
from analysis_driver import rest_communication
from analysis_driver.notification import default as ntf
from analysis_driver.exceptions import AnalysisDriverError, RestCommunicationError
from analysis_driver.app_logging import AppLogger
from analysis_driver.clarity import get_expected_yield_for_sample
from analysis_driver.constants import DATASET_NEW, DATASET_READY, DATASET_FORCE_READY, DATASET_PROCESSING,\
    DATASET_PROCESSED_SUCCESS, DATASET_PROCESSED_FAIL, DATASET_ABORTED, DATASET_REPROCESS, ELEMENT_RUN_NAME,\
    ELEMENT_NB_Q30_R2_CLEANED, ELEMENT_NB_Q30_R1_CLEANED, DATASET_DELETED

STATUS_VISIBLE = [DATASET_NEW, DATASET_REPROCESS, DATASET_READY, DATASET_FORCE_READY, DATASET_PROCESSING]
STATUS_HIDDEN = [DATASET_PROCESSED_SUCCESS, DATASET_PROCESSED_FAIL, DATASET_ABORTED, DATASET_DELETED]


class MostRecentProc:
    def __init__(self, dataset_type, dataset_name, initial_content=None):
        self.dataset_type = dataset_type
        self.dataset_name = dataset_name
        if initial_content:
            self._rest_entity = dict(initial_content)
            self.local_entity = dict(initial_content)
        else:
            self._rest_entity = None
            self.local_entity = dict(self.rest_entity)

    @property
    def rest_entity(self):
        if self._rest_entity is None:
            procs = rest_communication.get_documents(
                'analysis_driver_procs',
                where={'dataset_type': self.dataset_type, 'dataset_name': self.dataset_name},
                sort='-_created'
            )
            if procs:
                self._rest_entity = procs[0]
            else:
                self.initialise_entity()
        return self._rest_entity

    def initialise_entity(self):
        proc_id = '_'.join((self.dataset_type, self.dataset_name, self._now()))
        entity = {
            'proc_id': proc_id,
            'dataset_type': self.dataset_type,
            'dataset_name': self.dataset_name
        }
        rest_communication.post_entry('analysis_driver_procs', entity)
        rest_communication.patch_entry(
            self.dataset_type + 's',
            {'analysis_driver_procs': [proc_id]},
            id_field=self.dataset_type + '_id',
            element_id=self.dataset_name,
            update_lists=['analysis_driver_procs']
        )
        self._rest_entity = dict(entity)
        self.local_entity = dict(entity)

    def sync(self):
        patch_content = {}
        for k, v in self.local_entity.items():
            if v != self.rest_entity.get(k):
                patch_content[k] = v

        if patch_content:
            patch_success = rest_communication.patch_entry(
                'analysis_driver_procs',
                patch_content,
                id_field='proc_id',
                element_id=self.local_entity['proc_id']
            )
            if not patch_success:
                raise RestCommunicationError('Sync failed: ' + str(patch_content))
            self._rest_entity = dict(self.local_entity)

    def update_entity(self, **kwargs):
        self.local_entity.update(kwargs)
        self.sync()

    def change_status(self, status):
        self.update_entity(status=status)

    def start(self):
        self.update_entity(status=DATASET_PROCESSING, pid=os.getpid())

    def finish(self, status):
        self.update_entity(status=status, pid=None, end_date=self._now())

    def start_stage(self, stage_name):
        stages = self.local_entity.get('stages', [])
        new_stage = {'date_started': self._now(), 'stage_name': stage_name}
        stages.append(new_stage)
        self.update_entity(stages=stages)

    def end_stage(self, stage_name, exit_status=0):
        stages = self.local_entity['stages']
        for s in stages:
            if s['stage_name'] == stage_name:
                s.update({'date_finished': self._now(), 'exit_status': exit_status})

        self.update_entity(stages=stages)

    def get(self, key, ret_default=None):
        return self.local_entity.get(key, ret_default)

    @staticmethod
    def _now():
        return datetime.utcnow().strftime('%d_%m_%Y_%H:%M:%S')


class Dataset:
    type = None
    endpoint = None
    id_field = None

    def __init__(self, name, most_recent_proc=None):
        self.name = name
        if most_recent_proc:
            self.most_recent_proc = most_recent_proc
        else:
            self.most_recent_proc = MostRecentProc(self.type, self.name)

    @property
    def dataset_status(self):
        db_proc_status = self.most_recent_proc.get('status')
        if db_proc_status in (DATASET_REPROCESS, None):
            if self._is_ready():
                return DATASET_READY
            else:
                return DATASET_NEW
        else:
            return db_proc_status

    @property
    def stages(self):
        return [s['stage_name'] for s in self.most_recent_proc.get('stages', []) if 'date_finished' not in s]

    def start(self):
        assert self.dataset_status in (DATASET_READY, DATASET_FORCE_READY, DATASET_NEW)
        self.most_recent_proc.initialise_entity()  # take a new entity
        self.most_recent_proc.start()
        ntf.start_pipeline()

    def succeed(self, quiet=False):
        assert self.dataset_status == DATASET_PROCESSING
        self.most_recent_proc.finish(DATASET_PROCESSED_SUCCESS)
        if not quiet:
            ntf.end_pipeline(0)

    def fail(self, exit_status):
        assert self.dataset_status == DATASET_PROCESSING
        ntf.end_pipeline(exit_status)
        self.most_recent_proc.finish(DATASET_PROCESSED_FAIL)

    def abort(self):
        self.most_recent_proc.finish(DATASET_ABORTED)

    def reset(self):
        self.most_recent_proc.change_status(DATASET_REPROCESS)

    def start_stage(self, stage_name):
        ntf.start_stage(stage_name)
        self.most_recent_proc.start_stage(stage_name)

    def end_stage(self, stage_name, exit_status=0):
        ntf.end_stage(stage_name, exit_status)
        self.most_recent_proc.end_stage(stage_name, exit_status)

    @property
    def _is_ready(self):
        raise NotImplementedError

    def __str__(self):
        out = [self.name]
        pid = self.most_recent_proc.get('pid')
        if pid:
            out.append('(%s)' % pid)
        if self.stages:
            out.append('-- %s' % (', '.join(self.stages)))
        return ' '.join(out)

    __repr__ = __str__

    def __lt__(self, other):
        return self.name < other.name


class NoCommunicationDataset(Dataset):
    """Dummy dataset that can be used in QC object but won't contact the API"""
    type = "Notype"

    def start_stage(self, stage_name):
        pass

    def end_stage(self, stage_name, exit_status=0):
        pass

    def _change_status(self, status, finish=True):
        pass

    def _is_ready(self):
        pass


class RunDataset(Dataset):
    type = 'run'
    endpoint = 'runs'
    id_field = 'run_id'

    def __init__(self, name, path, use_int_dir, most_recent_proc=None):
        super().__init__(name, most_recent_proc)
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

    def __init__(self, name, most_recent_proc=None):
        super().__init__(name, most_recent_proc)
        self.run_elements = self._read_data()
        self._data_threshold = None

    def force(self):
        self.most_recent_proc.change_status(DATASET_FORCE_READY)

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


class DatasetScanner(AppLogger):
    endpoint = None
    item_id = None

    def __init__(self, config):
        self.input_dir = config.get('input_dir')
        self.__triggerignore = None

    def scan_datasets(self):
        return self._datasets_by_status(*STATUS_VISIBLE)

    def scan_processable_datasets(self):
        return self._datasets_by_status(
            DATASET_NEW,
            DATASET_REPROCESS,
            DATASET_READY,
            DATASET_FORCE_READY
        )

    def scan_hidden_datasets(self):
        return self._datasets_by_status(*STATUS_HIDDEN)

    def get_dataset(self, name, most_recent_proc=None):
        raise NotImplementedError

    def report(self, all_datasets=False):
        out = [
            '========= %s report =========' % self.__class__.__name__,
            'dataset location: ' + self.input_dir
        ]
        visible_datasets = self.scan_datasets()
        for k in sorted(visible_datasets):
            datasets = [str(d) for d in visible_datasets[k]]
            if datasets:
                out.append('=== ' + k + ' ===')
                out.append('\n'.join(datasets))

        if all_datasets:
            hidden_datasets = self.scan_hidden_datasets()
            for k in sorted(hidden_datasets):
                datasets = [str(d) for d in hidden_datasets[k]]
                if datasets:
                    out.append('=== ' + k + ' ===')
                    out.append('\n'.join(datasets))

        out.append('_' * 42)
        print('\n'.join(out))

    def _get_dataset_records_for_status(self, status):
        return [
            d for d in rest_communication.get_documents(
                self.endpoint,
                match={'proc_status': status}
            )
            if d[self.item_id] not in self._triggerignore
        ]

    def _get_datasets_for_status(self, status):
        return [
            self.get_dataset(d[self.item_id], d.get('most_recent_proc'))
            for d in self._get_dataset_records_for_status(status)
        ]

    def _datasets_by_status(self, *rest_api_statuses):
        datasets = defaultdict(list)
        for s in rest_api_statuses:
            for d in self._get_datasets_for_status(s):
                datasets[d.dataset_status].append(d)

        for k in datasets:
            datasets[k].sort()
        return datasets

    @property
    def _triggerignore(self):
        if self.__triggerignore is None:
            triggerignore = os.path.join(self.input_dir, '.triggerignore')

            ignorables = []
            if os.path.isfile(triggerignore):
                with open(triggerignore, 'r') as f:
                    for p in f.readlines():
                        if not p.startswith('#'):
                            ignorables.append(p.rstrip('\n'))
            self.debug('Ignoring %s datasets', len(ignorables))
            self.__triggerignore = ignorables
        return self.__triggerignore


class RunScanner(DatasetScanner):
    endpoint = 'aggregate/all_runs'
    item_id = 'run_id'
    expected_bcl_subdirs = ('SampleSheet.csv', 'RunInfo.xml', 'Data')

    def __init__(self, config):
        super().__init__(config)
        self.use_int_dir = 'intermediate_dir' in config

    def get_dataset(self, name, most_recent_proc=None):
        dataset_path = os.path.join(self.input_dir, name)
        if os.path.exists(dataset_path):
            return RunDataset(
                name,
                os.path.join(self.input_dir, name),
                use_int_dir=self.use_int_dir,
                most_recent_proc=most_recent_proc
            )

    def _get_dataset_records_for_status(self, status):
        if status == DATASET_NEW:
            datasets = rest_communication.get_documents(self.endpoint)
            for d in self._datasets_on_disk():
                if d not in datasets:
                    datasets.append({self.item_id: d})
            return datasets
        else:
            return super()._get_dataset_records_for_status(status)

    def _datasets_on_disk(self):
        return [
            d for d in os.listdir(self.input_dir)
            if self._is_valid_dataset(d) and d not in self._triggerignore
        ]

    def _is_valid_dataset(self, dataset):
        d = os.path.join(self.input_dir, dataset)
        if os.path.isdir(d) and not d.startswith('.'):
            observed = os.listdir(d)
            return all([subdir in observed for subdir in self.expected_bcl_subdirs])
        else:
            return False


class SampleScanner(DatasetScanner):
    endpoint = 'aggregate/samples'
    item_id = 'sample_id'

    def __init__(self, config):
        super().__init__(config)
        self.data_threshold = config.get('data_threshold')

    def get_dataset(self, name, most_recent_proc=None):
        return SampleDataset(name, most_recent_proc)
