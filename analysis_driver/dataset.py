import os
import re
from multiprocessing import Lock
from datetime import datetime
from errno import ESRCH
from os.path import join
from sys import modules
from time import sleep

from egcg_core import rest_communication, clarity
from egcg_core.app_logging import AppLogger
from egcg_core.config import cfg
from egcg_core.constants import *  # pylint: disable=unused-import
from egcg_core.exceptions import RestCommunicationError

from analysis_driver import reader
from analysis_driver.exceptions import AnalysisDriverError
from analysis_driver.notification import NotificationCentre
from analysis_driver.util import generate_samplesheet


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
            return DATASET_READY if self._is_ready() else DATASET_NEW
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
        self._assert_status(DATASET_READY, DATASET_FORCE_READY, DATASET_NEW, DATASET_RESUME)
        if self.dataset_status != DATASET_RESUME:
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

    def resume(self):
        self.terminate()
        self.most_recent_proc.change_status(DATASET_RESUME)

    def reset(self):
        self.terminate()
        self.most_recent_proc.change_status(DATASET_REPROCESS)

    def skip(self):
        self.most_recent_proc.finish(DATASET_PROCESSED_SUCCESS)

    def force(self):
        self.most_recent_proc.change_status(DATASET_FORCE_READY)

    def terminate(self):
        pid = self.most_recent_proc.get('pid')
        self.debug('Attempting to terminate pid %s for %s %s', pid, self.type, self.name)
        if not pid or not self._pid_valid(pid):
            self.debug('Attempted to terminate invalid pid %s', pid)
            return

        os.kill(pid, 10)
        while self._pid_running(pid):
            sleep(1)
        self.info('Terminated pid %s for %s %s', pid, self.type, self.name)

    def start_stage(self, stage_name):
        self.ntf.start_stage(stage_name)
        self.most_recent_proc.start_stage(stage_name)

    def end_stage(self, stage_name, exit_status=0):
        self.ntf.end_stage(stage_name, exit_status)
        self.most_recent_proc.end_stage(stage_name, exit_status)

    @staticmethod
    def _pid_valid(pid):
        cmd_file = os.path.join('/', 'proc', str(pid), 'cmdline')
        if os.path.isfile(cmd_file):
            with open(cmd_file, 'r') as f:
                return modules['__main__'].__file__ in f.read()
        return False

    @staticmethod
    def _pid_running(pid):
        try:
            os.kill(pid, 0)
        except OSError as err:
            if err.errno == ESRCH:  # no such process
                return False
        return True

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

    def __init__(self, name, most_recent_proc=None):
        super().__init__(name, most_recent_proc)
        self._run_info = None
        self._sample_sheet_file = None
        self.input_dir = os.path.join(cfg['input_dir'], self.name)
        self._run_elements = None
        self._barcode_len = None
        self._lims_run = None

    @property
    def run_info(self):
        if self._run_info is None:
            self._run_info = reader.RunInfo(self.input_dir)
        return self._run_info

    @property
    def sample_sheet_file(self):
        if self._sample_sheet_file is None:
            self._sample_sheet_file = join(self.input_dir, 'SampleSheet_analysis_driver.csv')
            generate_samplesheet(
                self,
                self._sample_sheet_file,
            )
        return self._sample_sheet_file

    @property
    def mask(self):
        return self.run_info.reads.generate_mask(self.barcode_len)

    def _is_ready(self):
        return True

    @property
    def run_elements(self):
        if not self._run_elements:
            self._run_elements = self._run_elements_from_lims()
        return self._run_elements

    def _run_elements_from_lims(self):
        run_elements = []

        def find_pooling_step_for_artifact(art, max_iteration=10, expected_pooling_step_name=None):
            nb_iteration = 0
            while len(art.input_artifact_list()) == 1:
                art = art.input_artifact_list()[0]
                if nb_iteration == max_iteration:
                    raise ValueError('Cannot find pooling step after %s iteraction' % max_iteration)
                nb_iteration += 1
            if expected_pooling_step_name and art.parent_process.type.name != expected_pooling_step_name:
                raise ValueError(
                    'Mismatching Step name: %s != %s' % (expected_pooling_step_name, art.parent_process.type.name)
                )
            return art.input_artifact_list()

        flowcell = set(self.lims_run.parent_processes()).pop().output_containers()[0]
        for lane in flowcell.placements:
            if len(flowcell.placements[lane].reagent_labels) > 1:
                artifacts = find_pooling_step_for_artifact(flowcell.placements[lane],
                                                           expected_pooling_step_name='Create PDP Pool')
            else:
                artifacts = [flowcell.placements[lane]]
            for artifact in artifacts:
                assert len(artifact.samples) == 1
                assert len(artifact.reagent_labels) == 1
                sample = artifact.samples[0]
                reagent_label = artifact.reagent_labels[0]
                match = re.match('(\w{4})-(\w{4}) \(([ATCG]{8})-([ATCG]{8})\)', reagent_label)
                run_element = {
                    ELEMENT_PROJECT_ID: sample.project.name,
                    ELEMENT_SAMPLE_INTERNAL_ID: sample.name,
                    ELEMENT_LIBRARY_INTERNAL_ID: sample.id,  # This is not the library id but it is unique
                    ELEMENT_LANE: lane.split(':')[0],
                    ELEMENT_BARCODE: ''
                }
                if self.has_barcodes:
                    run_element[ELEMENT_BARCODE] = match.group(3)
                run_elements.append(run_element)
        return run_elements

    @property
    def has_barcodes(self):
        return self.run_info.reads.has_barcodes

    @property
    def barcode_len(self):
        if not self._barcode_len:
            self._barcode_len = self._check_barcodes()
        return self._barcode_len

    def _check_barcodes(self):
        """
        For each run element, check that all the DNA barcodes are the same length
        :return: The DNA barcode length
        """
        previous_r = None

        for r in self.run_elements:
            if previous_r and len(previous_r[ELEMENT_BARCODE]) != len(r[ELEMENT_BARCODE]):
                raise AnalysisDriverError(
                    'Unexpected barcode length for %s: %s in project %s' % (
                        r[ELEMENT_SAMPLE_INTERNAL_ID], r[ELEMENT_BARCODE], r[ELEMENT_PROJECT_ID]
                    )
                )
            previous_r = r

        self.debug('Barcode check done. Barcode len: %s', len(r[ELEMENT_BARCODE]))
        return len(r[ELEMENT_BARCODE])

    @property
    def lims_run(self):
        if not self._lims_run:
            self._lims_run = clarity.get_run(self.name)
        return self._lims_run

    def is_sequencing(self):
        # assume that not being in the LIMS counts as 'is_sequencing'
        if not self.lims_run:
            self.warning('Run %s not found in the LIMS', self.name)
            return True
        # force the LIMS to update the RunStatus rather than passing the same cached RunStatus
        self.lims_run.get(force=True)
        return self.lims_run.udf.get('Run Status') in ['RunStarted', 'RunPaused']

    @property
    def run_metrics(self):
        return rest_communication.get_document('aggregate/all_runs', match={'run_id': self.name})

    @property
    def lane_metrics(self):
        return rest_communication.get_documents('aggregate/run_elements_by_lane', match={'run_id': self.name})


class SampleDataset(Dataset):
    type = 'sample'
    endpoint = 'samples'
    id_field = 'sample_id'

    def __init__(self, name, most_recent_proc=None, data_threshold=None):
        super().__init__(name, most_recent_proc)
        self._run_elements = None
        self._non_useable_run_elements = None
        self._data_threshold = data_threshold
        self._species = None
        self._genome_version = None
        self._user_sample_id = None

    @property
    def species(self):
        if self._species is None:
            self._species = clarity.get_species_from_sample(self.name)
        return self._species

    @property
    def genome_version(self):
        if self._genome_version is None:
            g = clarity.get_sample(self.name).udf.get('Genome Version')
            if g is None:
                g = cfg.query('species', self.species, 'default')
            self._genome_version = g
        return self._genome_version

    @property
    def reference_genome(self):
        return cfg['genomes'][self.genome_version]['fasta']

    @property
    def user_sample_id(self):
        if self._user_sample_id is None:
            self._user_sample_id = clarity.get_user_sample_name(self.name)
        return self._user_sample_id

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
        self._species = None
        self._genome_version = None

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

    @property
    def species(self):
        if self._species is None:
            s = set()
            for sample in self.samples_processed:
                species = sample.get('species_name', clarity.get_species_from_sample(sample['sample_id']))
                s.add(species)
            if len(s) != 1:
                raise AnalysisDriverError('Unexpected number of species (%s) in this project' % ', '.join(s))
            self._species = s.pop()
        return self._species

    @property
    def genome_version(self):
        if self._genome_version is None:
            g = clarity.get_sample(self.samples_processed[0]['sample_id']).udf.get('Genome Version')
            if not g:
                g = cfg.query('species', self.species, 'default')
            self._genome_version = g
        return self._genome_version

    @property
    def reference_genome(self):
        return cfg['genomes'][self.genome_version]['fasta']

    def __str__(self):
        return '%s  (%s samples / %s) ' % (super().__str__(), len(self.samples_processed), self.number_of_samples)


class MostRecentProc:
    def __init__(self, dataset_type, dataset_name, initial_content=None):
        self.lock = Lock()
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
        doc = rest_communication.get_document(
            'analysis_driver_stages',
            where={'stage_id': self._stage_id(stage_name)}
        )
        if doc:
            rest_communication.patch_entry(
                'analysis_driver_stages',
                {'date_started': self._now(), 'date_finished': None, 'exit_status': None},
                'stage_id', self._stage_id(stage_name)
            )
        else:
            rest_communication.post_entry(
                'analysis_driver_stages',
                {'stage_id': self._stage_id(stage_name), 'date_started': self._now(),
                 'stage_name': stage_name, 'analysis_driver_proc': self.proc_id}
            )
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
