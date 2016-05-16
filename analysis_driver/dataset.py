import os
import threading
from datetime import datetime
from analysis_driver import rest_communication
from analysis_driver.notification import default as ntf
from analysis_driver.clarity import get_expected_yield_for_sample
from analysis_driver.exceptions import AnalysisDriverError, RestCommunicationError
from analysis_driver.constants import DATASET_NEW, DATASET_READY, DATASET_FORCE_READY, DATASET_REPROCESS,\
    DATASET_PROCESSING, DATASET_PROCESSED_SUCCESS, DATASET_PROCESSED_FAIL, DATASET_ABORTED, ELEMENT_RUN_NAME,\
    ELEMENT_NB_Q30_R1_CLEANED, ELEMENT_NB_Q30_R2_CLEANED


class Dataset:
    type = None
    endpoint = None
    id_field = None

    def __init__(self, name, most_recent_proc=None):
        self.name = name
        self.most_recent_proc = MostRecentProc(self.type, self.name, most_recent_proc)

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
        self._assert_status(DATASET_READY, DATASET_FORCE_READY, DATASET_NEW, method='start')
        self.most_recent_proc.initialise_entity()  # take a new entity
        self.most_recent_proc.start()
        ntf.start_pipeline()

    def succeed(self, quiet=False):
        self._assert_status(DATASET_PROCESSING, method='succeed')
        self.most_recent_proc.finish(DATASET_PROCESSED_SUCCESS)
        if not quiet:
            ntf.end_pipeline(0)

    def fail(self, exit_status):
        self._assert_status(DATASET_PROCESSING, method='fail')
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

    def _assert_status(self, *allowed_statuses, method=None):
        status = self.dataset_status
        assert status in allowed_statuses, 'Tried to %s a %s %s' % (method, status, self.type)

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

    def __init__(self, name, most_recent_proc=None, data_threshold=None):
        super().__init__(name, most_recent_proc)
        self._run_elements = None
        self._data_threshold = data_threshold

    def force(self):
        self.most_recent_proc.change_status(DATASET_FORCE_READY)

    @property
    def run_elements(self):
        if self._run_elements is None:
            self._run_elements = self._read_data()
        return self._run_elements

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


class MostRecentProc:
    def __init__(self, dataset_type, dataset_name, initial_content=None):
        self.lock = threading.Lock()
        self.dataset_type = dataset_type
        self.dataset_name = dataset_name
        if initial_content:
            self._entity = initial_content.copy()
        else:
            self._entity = None

    @property
    def entity(self):
        if self._entity is None:
            procs = rest_communication.get_documents(
                'analysis_driver_procs',
                where={'dataset_type': self.dataset_type, 'dataset_name': self.dataset_name},
                sort='-_created'
            )
            if procs:
                #Remove the private (starting with _ ) fields from the dict
                self._entity = {k: v for k, v in procs[0].items() if not k.startswith('_')}
            else:
                self.initialise_entity()
        return self._entity

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
        self._entity = entity

    def sync(self):
        patch_content = self.entity.copy()
        element_id = patch_content.pop('proc_id')
        patch_success = rest_communication.patch_entry(
            'analysis_driver_procs',
            patch_content,
            id_field='proc_id',
            element_id=element_id
        )
        if not patch_success:
            raise RestCommunicationError('Sync failed: ' + str(patch_content))

    def update_entity(self, **kwargs):
        with self.lock:
            self.entity.update(kwargs)
            self.sync()

    def change_status(self, status):
        self.update_entity(status=status)

    def start(self):
        self.update_entity(status=DATASET_PROCESSING, pid=os.getpid())

    def finish(self, status):
        self.update_entity(status=status, pid=None, end_date=self._now())

    def start_stage(self, stage_name):
        stages = self.entity.get('stages', [])
        new_stage = {'date_started': self._now(), 'stage_name': stage_name}
        stages.append(new_stage)
        self.update_entity(stages=stages)

    def end_stage(self, stage_name, exit_status=0):
        stages = self.entity['stages']
        for s in stages:
            if s['stage_name'] == stage_name:
                s.update({'date_finished': self._now(), 'exit_status': exit_status})

        self.update_entity(stages=stages)

    def get(self, key, ret_default=None):
        return self.entity.get(key, ret_default)

    @staticmethod
    def _now():
        return datetime.utcnow().strftime('%d_%m_%Y_%H:%M:%S')
