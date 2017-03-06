import os
import threading
from sys import modules
from datetime import datetime
from analysis_driver.notification import NotificationCentre
from egcg_core.app_logging import logging_default as log_cfg
from egcg_core import rest_communication, clarity
from egcg_core.app_logging import AppLogger
from egcg_core.exceptions import RestCommunicationError
from analysis_driver.exceptions import AnalysisDriverError
from egcg_core.constants import DATASET_NEW, DATASET_READY, DATASET_FORCE_READY, DATASET_REPROCESS,\
    DATASET_PROCESSING, DATASET_PROCESSED_SUCCESS, DATASET_PROCESSED_FAIL, DATASET_ABORTED, ELEMENT_RUN_NAME,\
    ELEMENT_NB_Q30_R1_CLEANED, ELEMENT_NB_Q30_R2_CLEANED

app_logger = log_cfg.get_logger('versions')

class Dataset(AppLogger):
    type = None
    endpoint = None
    id_field = None

    def __init__(self, name, most_recent_proc=None):
        self.name = name
        self.most_recent_proc = MostRecentProc(self.type, self.name, most_recent_proc)
        self._ntf = None

    @property
    def ntf(self):
        if self._ntf is None:
            self._ntf = NotificationCentre(self.name)
        return self._ntf

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
    def running_stages(self):
        stages = rest_communication.get_documents(
            'analysis_driver_stages',
            all_pages=True,
            where={'analysis_driver_proc': self.most_recent_proc.get('proc_id'), 'date_finished': None}
        )
        return [s['stage_name'] for s in stages]

    def get_stage(self, stage_name):
        return rest_communication.get_document(
            'analysis_driver_stages',
            where={'analysis_driver_proc': self.most_recent_proc.get('proc_id'), 'stage_name': stage_name}
        )

    def start(self):
        self._assert_status(DATASET_READY, DATASET_FORCE_READY, DATASET_NEW)
        self.most_recent_proc.initialise_entity()  # take a new entity
        self.most_recent_proc.start()
        self.ntf.start_pipeline()

    def succeed(self):
        self._assert_status(DATASET_PROCESSING)
        self.most_recent_proc.finish(DATASET_PROCESSED_SUCCESS)
        self.ntf.end_pipeline(0)

    def fail(self, exit_status):
        self._assert_status(DATASET_PROCESSING)
        self.ntf.end_pipeline(exit_status)
        self.most_recent_proc.finish(DATASET_PROCESSED_FAIL)

    def abort(self):
        self.most_recent_proc.finish(DATASET_ABORTED)

    def reset(self):
        self.terminate()
        self.most_recent_proc.change_status(DATASET_REPROCESS)

    def skip(self):
        self.most_recent_proc.finish(DATASET_PROCESSED_SUCCESS)

    def force(self):
        self.most_recent_proc.change_status(DATASET_FORCE_READY)

    def terminate(self):
        pid = self.most_recent_proc.get('pid')
        if pid and self._is_valid_pid(pid):
            self.info('Terminating pid %s for %s %s', pid, self.type, self.name)
            os.kill(pid, 10)

    def start_stage(self, stage_name):
        self.ntf.start_stage(stage_name)
        self.most_recent_proc.start_stage(stage_name)

    def end_stage(self, stage_name, exit_status=0):
        self.ntf.end_stage(stage_name, exit_status)
        self.most_recent_proc.end_stage(stage_name, exit_status)

    @staticmethod
    def _is_valid_pid(pid):
        cmd_file = os.path.join('/', 'proc', str(pid), 'cmdline')
        if os.path.isfile(cmd_file):
            with open(cmd_file, 'r') as f:
                return modules['__main__'].__file__ in f.read()

    def _assert_status(self, *allowed_statuses):
        # make sure the most recent process is the same as the one in the REST API
        self.most_recent_proc.retrieve_entity()
        status = self.dataset_status
        assert status in allowed_statuses, 'Status assertion failed on a %s %s' % (status, self.type)

    @property
    def _is_ready(self):
        raise NotImplementedError

    def __str__(self):
        s = self.name
        pid = self.most_recent_proc.get('pid')
        if pid:
            s += ' (%s)' % pid
        stages = self.running_stages
        if stages:
            s += ' -- ' + ', '.join(stages)
        return s

    __repr__ = __str__

    def __lt__(self, other):
        return self.name < other.name


class NoCommunicationDataset(Dataset):
    """Dummy dataset that can be used in QC object but won't contact the API"""
    type = 'Notype'

    def start_stage(self, stage_name):
        pass

    def end_stage(self, stage_name, exit_status=0):
        pass

    def _change_status(self, status, finish=True):
        pass

    def _is_ready(self):
        pass

    def __str__(self):
        return self.name

    __repr__ = __str__


class RunDataset(Dataset):
    type = 'run'
    endpoint = 'runs'
    id_field = 'run_id'

    def __init__(self, name, path, most_recent_proc=None):
        super().__init__(name, most_recent_proc)
        self.path = path

    def _is_ready(self):
        return True

    def is_sequencing(self):
        # Assume the run has started and not finished if the status is 'RunStarted' or if it hasn't yet appeared in the LIMS
        if not clarity.get_run(self.name):
            app_logger.warning('Run %s not found in the LIMS', self.name)
            return True
        return clarity.get_run(self.name).udf.get('Run Status') == 'RunStarted'


class SampleDataset(Dataset):
    type = 'sample'
    endpoint = 'samples'
    id_field = 'sample_id'

    def __init__(self, name, most_recent_proc=None, data_threshold=None):
        super().__init__(name, most_recent_proc)
        self._run_elements = None
        self._non_useable_run_elements = None
        self._data_threshold = data_threshold

    @property
    def run_elements(self):
        if self._run_elements is None:
            self._run_elements = rest_communication.get_documents(
                'run_elements', where={'sample_id': self.name, 'useable': 'yes'}
            )
        return self._run_elements

    @property
    def non_useable_run_elements(self):
        if self._non_useable_run_elements is None:
            self._non_useable_run_elements = rest_communication.get_documents(
                'run_elements', where={'sample_id': self.name, 'useable': {'$ne': 'yes'}}
            )
        return self._non_useable_run_elements

    def _amount_data(self):
        return sum(
            [
                int(r.get(ELEMENT_NB_Q30_R1_CLEANED, 0)) + int(r.get(ELEMENT_NB_Q30_R2_CLEANED, 0))
                for r in self.run_elements
            ]
        )

    @property
    def data_threshold(self):
        if self._data_threshold is None:
            self._data_threshold = clarity.get_expected_yield_for_sample(self.name)
        if not self._data_threshold:
            raise AnalysisDriverError('Could not find data threshold in LIMS for ' + self.name)
        return self._data_threshold

    def _is_ready(self):
        return self.data_threshold and int(self._amount_data()) > int(self.data_threshold)

    def __str__(self):
        runs = sorted(set(r.get(ELEMENT_RUN_NAME) for r in self.run_elements))
        non_useable_runs = sorted(set(r.get(ELEMENT_RUN_NAME) for r in self.non_useable_run_elements))

        s = '%s  (%s / %s  from %s) ' % (
            super().__str__(), self._amount_data(), self.data_threshold, ', '.join(runs)
        )
        if non_useable_runs:
            s += '(non useable run elements in %s)' % ', '.join(non_useable_runs)
        else:
            s += '(no non useable run elements)'
        return s


class ProjectDataset(Dataset):
    type = 'project'
    endpoint = 'projects'
    id_field = 'project_id'

    def __init__(self, name, most_recent_proc=None):
        super().__init__(name, most_recent_proc)
        self._number_of_samples = None
        self._samples_processed = None

    def _is_ready(self):
        if self.number_of_samples > 0 and self.number_of_samples > len(self.samples_processed):
            return True
        else:
            return False

    @property
    def samples_processed(self):
        if not self._samples_processed:
            self._samples_processed = rest_communication.get_documents(
                'aggregate/samples',
                match={'project_id': self.name, 'proc_status': 'finished'}
            )
        return self._samples_processed

    @property
    def number_of_samples(self):
        if not self._number_of_samples:
            project_from_lims = clarity.get_project(self.name)
            if project_from_lims:
                self._number_of_samples = project_from_lims.udf.get('Number of Quoted Samples')
                if not self._number_of_samples:
                    self._number_of_samples = -1
            else:
                raise AnalysisDriverError('Could not find number of quoted samples in LIMS for ' + self.name)
        return self._number_of_samples

    def __str__(self):
        return '%s  (%s samples / %s) ' % ( super().__str__(), len(self.samples_processed), self.number_of_samples)


class MostRecentProc:
    def __init__(self, dataset_type, dataset_name, initial_content=None):
        self.lock = threading.Lock()
        self.dataset_type = dataset_type
        self.dataset_name = dataset_name
        self.proc_id = None
        if initial_content:
            self._entity = initial_content.copy()
        else:
            self._entity = None

    def retrieve_entity(self):
        procs = rest_communication.get_documents(
            'analysis_driver_procs',
            where={'dataset_type': self.dataset_type, 'dataset_name': self.dataset_name},
            sort='-_created'
        )
        if procs:
            # remove the private (starting with '_') fields from the dict
            self._entity = {k: v for k, v in procs[0].items() if not k.startswith('_')}
        else:
            self._entity = {}

    @property
    def entity(self):
        if self._entity is None:
            self.retrieve_entity()
        return self._entity

    def initialise_entity(self):
        self.proc_id = '_'.join((self.dataset_type, self.dataset_name, self._now()))
        entity = {
            'proc_id': self.proc_id,
            'dataset_type': self.dataset_type,
            'dataset_name': self.dataset_name
        }
        rest_communication.post_entry('analysis_driver_procs', entity)
        rest_communication.patch_entry(
            self.dataset_type + 's',
            {'analysis_driver_procs': [self.proc_id]},
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
            if not self.entity:  # initialise self._entity with database content, or {} if none available
                self.initialise_entity()  # if self._entity == {}, then initialise and push to database

            self.entity.update(kwargs)
            self.sync()

    def change_status(self, status):
        self.update_entity(status=status)

    def start(self):
        self.update_entity(status=DATASET_PROCESSING, pid=os.getpid())

    def finish(self, status):
        self.update_entity(status=status, pid=None, end_date=self._now())

    def start_stage(self, stage_name):
        success = rest_communication.post_entry(
            'analysis_driver_stages',
            {'stage_id': self._stage_id(stage_name), 'date_started': self._now(), 'stage_name': stage_name,
             'analysis_driver_proc': self.proc_id}
        )
        assert success

        stages = self.entity.get('stages', [])
        stages.append(self._stage_id(stage_name))
        self.update_entity(stages=stages)

    def end_stage(self, stage_name, exit_status=0):
        rest_communication.patch_entry(
            'analysis_driver_stages',
            {'date_finished': self._now(), 'exit_status': exit_status}, 'stage_id', self._stage_id(stage_name)
        )

    def get(self, key, ret_default=None):
        return self.entity.get(key, ret_default)

    def _stage_id(self, stage_name):
        return self.proc_id + '_' + stage_name

    @staticmethod
    def _now():
        return datetime.utcnow().strftime('%d_%m_%Y_%H:%M:%S')
