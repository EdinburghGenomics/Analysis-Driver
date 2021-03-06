import os
import re
import signal
from multiprocessing import Lock, Queue
from datetime import datetime
from errno import ESRCH
from sys import modules
from time import sleep
from egcg_core import rest_communication, clarity
from egcg_core.app_logging import AppLogger
from egcg_core.clarity import sanitize_user_id, get_species_name
from egcg_core.config import cfg
from egcg_core.util import query_dict, find_file
from egcg_core.constants import *  # pylint: disable=unused-import
from egcg_core.exceptions import RestCommunicationError
from analysis_driver import reader
from analysis_driver.exceptions import AnalysisDriverError
from analysis_driver.notification import NotificationCentre, LimsNotification
from analysis_driver.tool_versioning import toolset
from analysis_driver.util.helper_functions import prepend_path_to_data_files
from analysis_driver.pipelines import demultiplexing, bcbio, qc, variant_calling, projects, variant_calling_gatk4,\
    qc_gatk4, human_variant_calling_gatk4


def now(datefmt='%d_%m_%Y_%H:%M:%S'):
    return datetime.now().strftime(datefmt)


class Dataset(AppLogger):
    type = None
    endpoint = None
    id_field = None
    exceptions = Queue()  # This needs to be a Queue to be accessible from every subprocess

    def __init__(self, name, most_recent_proc=None):
        self.name = name
        self.most_recent_proc = MostRecentProc(self, most_recent_proc)
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
    def data_source(self):
        return None

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

    def initialise_entity(self):
        self.most_recent_proc.initialise_entity()

    def start(self):
        self._assert_status(DATASET_READY, DATASET_FORCE_READY, DATASET_NEW, DATASET_RESUME)
        if self.dataset_status == DATASET_RESUME:
            self.most_recent_proc.retrieve_entity()  # use existing entity
        else:
            self.initialise_entity()  # take a new entity
        self.most_recent_proc.start()
        self.ntf.start_pipeline()

    def succeed(self):
        self._assert_status(DATASET_PROCESSING)
        self.most_recent_proc.finish(DATASET_PROCESSED_SUCCESS)
        self.ntf.end_pipeline(0)

    def fail(self, exit_status):
        self._assert_status(DATASET_PROCESSING)
        self.most_recent_proc.finish(DATASET_PROCESSED_FAIL)
        self.ntf.end_pipeline(exit_status)

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
        self._terminate(signal.SIGUSR2)

    def soft_terminate(self):
        """Send SIGUSR1 to analysis driver, causing Luigi to stop submitting any new jobs."""
        self._terminate(signal.SIGUSR1)

    def _terminate(self, signal_id):
        pid = self.most_recent_proc.get('pid')
        self.info('Attempting to terminate pid %s for %s %s with signal %s', pid, self.type, self.name, signal_id)
        if not pid or not self._pid_valid(pid):
            self.warning('Termination unsuccessful')
            return

        os.kill(pid, signal_id)
        while self._pid_running(pid):
            sleep(1)

        self.info('Termination successful')

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

    def register_exception(self, luigi_task, exception):
        self.exceptions.put((luigi_task.stage_name, exception))

    def raise_exceptions(self):
        if not self.exceptions.empty():
            first_exception = None
            while not self.exceptions.empty():
                name, exception = self.exceptions.get()
                if not first_exception:
                    first_exception = exception
                self.critical(
                    'exception: %s%s raised in %s',
                    exception.__class__.__name__,
                    exception.args,
                    name
                )
            # Only raise the first exception raised by tasks
            raise first_exception

    def __str__(self):
        return '%s(name=%s)' % (self.__class__.__name__, self.name)

    def report(self):
        s = self.name
        pid = self.most_recent_proc.get('pid')
        if pid:
            s += ' (%s)' % pid
        stages = self.running_stages
        if stages:
            s += ' -- ' + ', '.join(stages)
        return s

    def __lt__(self, other):
        return self.name < other.name

    @property
    def pipeline(self):
        return None

    def _processing_instruction(self):
        pipeline = self.pipeline
        return {
            'name': pipeline.name,
            'toolset_type': pipeline.toolset_type,
            'toolset_version': toolset.latest_version(pipeline.toolset_type)
        }

    def resolve_pipeline_and_toolset(self):
        instruction = self._processing_instruction()
        toolset.configure(
            instruction['toolset_type'],
            instruction['toolset_version'],
            os.path.join(cfg['jobs_dir'], self.name, 'program_versions.yaml')
        )


class NoCommunicationDataset(Dataset):
    """Dummy dataset that can be used in QC object but won't contact the API"""
    type = 'Notype'

    def start_stage(self, stage_name):
        pass

    def end_stage(self, stage_name, exit_status=0):
        pass

    def _is_ready(self):
        pass


class NoCommunicationSampleDataset(NoCommunicationDataset):
    def __init__(self, name):
        super().__init__(name)
        self._run_elements = None

    @property
    def run_elements(self):
        if self._run_elements is None:
            self._run_elements = rest_communication.get_documents(
                'run_elements', where={'sample_id': self.name, 'useable': 'yes'}
            )
        return self._run_elements

    @property
    def project_id(self):
        return self.run_elements[0]['project_id']

    def _processing_instruction(self):
        instruction = rest_communication.get_document(
            'projects', where={'project_id': self.project_id}
        ).get('sample_pipeline')

        if not instruction:
            raise AnalysisDriverError('No instruction set in project %s' % self.project_id)

        return instruction


class RunDataset(Dataset):
    type = 'run'
    endpoint = 'runs'
    id_field = 'run_id'

    def __init__(self, name, most_recent_proc=None):
        super().__init__(name, most_recent_proc)
        self._run_info = None
        self._sample_sheet_file_per_lane = {}
        self.input_dir = os.path.join(cfg['run']['input_dir'], self.name)
        self._run_status = None
        self._run_elements = None
        self._barcode_len = {}
        self._rapid_samples_by_lane = None

    def initialise_entity(self):
        run = rest_communication.get_document('runs', where={'run_id': self.name})
        if not run:
            rest_communication.post_entry('runs', {'run_id': self.name})

        super().initialise_entity()

    @property
    def run_info(self):
        if self._run_info is None:
            self._run_info = reader.RunInfo(self.input_dir)
        return self._run_info

    def sample_sheet_file(self):
        pass

    def sample_sheet_file_for_lane(self, lane):
        if lane not in self._sample_sheet_file_per_lane:
            self._sample_sheet_file_per_lane[lane] = os.path.join(self.input_dir, 'SampleSheet_analysis_driver_lane%s.csv' % lane)
            if not os.path.isfile(self._sample_sheet_file_per_lane[lane]):
                self._generate_samplesheet(self._sample_sheet_file_per_lane[lane], lane)
        return self._sample_sheet_file_per_lane[lane]

    def _generate_samplesheet(self, filename, lane=None):
        all_lines = [
            '[Header]', 'Date, ' + now('%d/%m/%Y'), 'Workflow, Generate FASTQ Only', '',
            '[Settings]', 'Adapter, AGATCGGAAGAGCACACGTCTGAACTCCAGTCA',
            'AdapterRead2, AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT', '', '[Data]',
            'Lane,Sample_ID,Sample_Name,Sample_Project,index'
        ]
        # order the run elements so they produce a deterministic samplesheet
        for run_element in sorted(self.run_elements, key=lambda x: (x[ELEMENT_LANE], x[ELEMENT_BARCODE])):
            if lane and run_element[ELEMENT_LANE] != str(lane):
                continue

            all_lines.append(','.join([
                run_element[ELEMENT_LANE],
                run_element[ELEMENT_SAMPLE_INTERNAL_ID],
                run_element[ELEMENT_LIBRARY_INTERNAL_ID],
                run_element[ELEMENT_PROJECT_ID],
                run_element[ELEMENT_BARCODE]
            ]))
        with open(filename, 'w') as f:
            f.write('\n'.join(all_lines) + '\n')

    def mask_per_lane(self, lane):
        if self.has_barcode_in_lane(lane):
            return self.run_info.reads.generate_mask(self.barcode_len(lane))
        return self.run_info.reads.generate_mask(0)

    def _is_ready(self):
        return True

    @property
    def run_elements(self):
        if not self._run_elements:
            self._run_elements = self._run_elements_from_lims_endpoint()
        return self._run_elements

    @staticmethod
    def _find_pooling_step_for_artifact(art, expected_pooling_step_names, max_iterations=10):
        n = 1
        while len(art.input_artifact_list()) == 1:
            art = art.input_artifact_list()[0]
            if n >= max_iterations:
                raise ValueError('Cannot find pooling step after %s iterations' % max_iterations)
            n += 1
        if art.parent_process.type.name not in expected_pooling_step_names:
            raise ValueError('Unexpected step name: %s' % art.parent_process.type.name)
        return art.input_artifact_list()

    @property
    def number_of_lane(self):
        return len(self.run_status['lanes'])

    @property
    def rapid_samples_by_lane(self):
        """
        Search the sequencing run in the Lims for any samples with the UDF 'Rapid Analysis' set, and return their sample
        ID, project ID and all UDFs.
        """
        if self._rapid_samples_by_lane is None:
            self._rapid_samples_by_lane = {}
            for lane in self.run_status['lanes']:
                if len(lane['samples']) > 1:
                    continue  # we don't want to run rapid processing on pools
                sample_data = lane['samples'][0]

                sample_info = rest_communication.get_document('lims/sample_info',
                                                              match={'sample_id': sample_data['sample_id']})
                if sample_info.get('Rapid Analysis') == 'Yes':
                    self._rapid_samples_by_lane[str(lane['lane'])] = sample_info

        return self._rapid_samples_by_lane

    @property
    def run_status(self):
        if self._run_status is None:
            self._run_status = rest_communication.get_document('lims/run_status', match={'run_id': self.name})
            if not self._run_status:
                raise AnalysisDriverError('Run status information for %s is not available in the LIMS' % self.name)
        return self._run_status

    def _run_elements_from_lims_endpoint(self):
        run_elements = []

        for lane in self.run_status['lanes']:
            for sample_data in lane['samples']:
                run_element = {
                    ELEMENT_PROJECT_ID: sample_data['project_id'],
                    ELEMENT_SAMPLE_INTERNAL_ID: sample_data['sample_id'],
                    ELEMENT_LIBRARY_INTERNAL_ID: sample_data['artifact_id'],  # This is not the library id but it is unique
                    ELEMENT_LANE: str(lane['lane']),
                    ELEMENT_BARCODE: ''
                }
                if len(lane['samples']) > 1:
                    barcode = ''
                    for pattern in (
                            # TruSeq label, e.g, A412-A208 (ATGCATGC-CTGACTGA)
                            '(\w{4})-(\w{4}) \(([ATCG]{8})-([ATCG]{8})\)',
                            # IDT label, e.g, 001A IDT-ILMN TruSeq DNA-RNA UD 96 Indexes  Plate_UDI0001 (ATGCATGC-CTGACTGA)
                            '(\w{4}) IDT-ILMN TruSeq DNA-RNA UD 96 Indexes\s+(Plate_\w{7}) \(([ATGC]{8})-([ATGC]{8})\)'
                    ):
                        match = re.match(pattern, sample_data['barcode'])
                        if match:
                            barcode = match.group(3)

                    if not barcode:
                        raise AnalysisDriverError('Invalid barcode found: %s' % sample_data['barcode'])
                    run_element[ELEMENT_BARCODE] = barcode
                run_elements.append(run_element)

        return run_elements

    def has_barcode_in_lane(self, lane_number):
        """
        Returns False if this lane has only one barcode (regardless of what is in the RunInfo.xml).
        Otherwise delegate to the RunInfo object.
        :param lane_number: the lane number to query.
        :return: True if this lane needs to be use the barcode False otherwise.
        """
        lane_data = [l for l in self.run_status['lanes'] if str(l['lane']) == str(lane_number)][0]
        # If there is only one sample in this lane then we won't demultiplex and should use a barcode length of 0
        if len(lane_data['samples']) == 1:
            return False
        return self.run_info.reads.has_barcodes

    def barcode_len(self, lane):
        if lane not in self._barcode_len:
            self._barcode_len[lane] = self._check_barcodes(lane)
        return self._barcode_len[lane]

    def _check_barcodes(self, lane):
        """
        For each run element of a lane, check that all the DNA barcodes are the same length
        :param lane: the lane for which the barcode len should be returned
        :return: The DNA barcode length for the specified lane
        """
        previous_r = None

        for r in self.run_elements:
            if r[ELEMENT_LANE] != str(lane):
                continue
            if previous_r and len(previous_r[ELEMENT_BARCODE]) != len(r[ELEMENT_BARCODE]):
                raise AnalysisDriverError(
                    'Unexpected barcode length for %s: %s in project %s' % (
                        r[ELEMENT_SAMPLE_INTERNAL_ID], r[ELEMENT_BARCODE], r[ELEMENT_PROJECT_ID]
                    )
                )
            previous_r = r

        self.debug('Barcode check for lane %s done. Barcode len: %s', lane, len(previous_r[ELEMENT_BARCODE]))
        return len(previous_r[ELEMENT_BARCODE])

    def is_sequencing(self):
        # assume that not being in the LIMS counts as 'is_sequencing'
        if not self.run_status:
            self.warning('Run %s not found in the LIMS', self.name)
            return True
        # force the update of the RunStatus rather than passing the same cached values
        self._run_status = None
        return self.run_status.get('run_status') in ['RunStarted', 'RunPaused']

    @property
    def lane_metrics(self):
        return rest_communication.get_documents('lanes', where={'run_id': self.name})

    @property
    def pipeline(self):
        return demultiplexing.Demultiplexing(self)

    @staticmethod
    def reference_genome(run_element):
        return SampleDataset(run_element['sample_id']).reference_genome


class SampleDataset(Dataset):
    type = 'sample'
    endpoint = 'samples'
    id_field = 'sample_id'

    def __init__(self, name, most_recent_proc=None):
        super().__init__(name, most_recent_proc)
        self._run_elements = None
        self._sample = None
        self._lims_sample_info = None
        self._lims_sample_status = None
        self._non_useable_run_elements = None
        self._species = None
        self._genome_version = None
        self._user_sample_id = None
        self._lims_ntf = None
        self._genome_dict = None

    @property
    def lims_sample_info(self):
        if self._lims_sample_info is None:
            self._lims_sample_info = rest_communication.get_document('lims/sample_info', match={'sample_id': self.name})
            if self._lims_sample_info is None:
                self._lims_sample_info = {}
        return self._lims_sample_info

    @property
    def lims_ntf(self):
        if self._lims_ntf is None:
            self._lims_ntf = LimsNotification(self.name)
        return self._lims_ntf

    @property
    def species(self):
        if self._species is None:
            self._species = self.sample.get('species_name')
        if self._species is None:
            self._species = clarity.get_species_from_sample(self.name)
        return self._species

    @property
    def genome_version(self):
        if self._genome_version is None:
            self._genome_version = clarity.get_genome_version(sample_id=self.name, species=self.species)
        return self._genome_version

    @property
    def genome_dict(self):
        if self._genome_dict is None:
            # Getting reference genome data from rest API
            genome_response = rest_communication.get_document('genomes', where={'assembly_name': self.genome_version})
            # Checking project whitelist to ensure reference genome can be used
            if 'project_whitelist' in genome_response and self.project_id not in genome_response['project_whitelist']:
                raise AnalysisDriverError('Project ID ' + self.project_id + ' not in whitelist for reference genome '
                                          + self.genome_version)
            # Appending genomes_dir to data_files items
            genome_response['data_files'] = prepend_path_to_data_files(
                cfg.get('genomes_dir', ''), genome_response['data_files']
            )
            self._genome_dict = genome_response
        return self._genome_dict

    @property
    def reference_genome(self):
        return self.genome_dict['data_files']['fasta']

    @property
    def data_source(self):
        return [r['run_element_id'] for r in self.run_elements]

    @property
    def user_sample_id(self):
        return clarity.get_user_sample_name(self.name)

    @property
    def run_elements(self):
        if self._run_elements is None:
            self._run_elements = rest_communication.get_documents(
                'run_elements', where={'sample_id': self.name, 'useable': 'yes'}
            )
        return self._run_elements

    @property
    def sample(self):
        if self._sample is None:
            self._sample = rest_communication.get_document('samples', where={'sample_id': self.name})
        return self._sample

    @property
    def lims_sample_status(self):
        if self._lims_sample_status is None:
            self._lims_sample_status = rest_communication.get_document('lims/sample_status',
                                                                       match={'sample_id': self.name})
        return self._lims_sample_status

    def library_preparation(self):
        return self.lims_sample_status.get('library_type')

    @property
    def non_useable_run_elements(self):
        if self._non_useable_run_elements is None:
            self._non_useable_run_elements = rest_communication.get_documents(
                'run_elements', where={'sample_id': self.name, 'useable': {'$ne': 'yes'}}
            )
        return self._non_useable_run_elements

    @property
    def project_id(self):
        return self.sample['project_id']

    def _amount_data(self):
        y = query_dict(self.sample, 'aggregated.clean_yield_in_gb')
        if y:
            return int(y * 1000000000)
        else:
            return 0

    def _amount_coverage(self):
        return query_dict(self.sample, 'aggregated.from_run_elements.mean_coverage') or 0

    @property
    def pc_q30(self):
        return self.sample.get('aggregated', {}).get('clean_pc_q30') or 0

    @property
    def required_yield_threshold(self):
        if 'required_yield' in self.sample:
            return self.sample.get('required_yield')
        raise AnalysisDriverError('Could not find required yield threshold for ' + self.name)

    @property
    def required_coverage_threshold(self):
        if 'required_coverage' in self.sample:
            return self.sample.get('required_coverage')
        raise AnalysisDriverError('Could not find required coverage threshold for ' + self.name)

    def _processing_instruction(self):
        instruction = rest_communication.get_document(
            'projects', where={'project_id': self.project_id}
        ).get('sample_pipeline')

        if not instruction:
            instruction = super()._processing_instruction()
            rest_communication.patch_entry('projects', {'sample_pipeline': instruction}, 'project_id', self.project_id)

        return instruction

    def _is_ready(self):
        return self.pc_q30 > 75 and (
            (self.required_yield_threshold and self._amount_data() > self.required_yield_threshold) or
            (self.required_coverage_threshold and self._amount_coverage() > self.required_coverage_threshold)
        )

    def report(self):
        runs = query_dict(self.sample, 'aggregated.run_ids')
        non_useable_runs = sorted(set(r[ELEMENT_RUN_NAME] for r in self.non_useable_run_elements))

        s = '%s  (yield: %s / %s and coverage: %s / %s from %s) ' % (
            super().report(), self._amount_data(), self.required_yield_threshold, self._amount_coverage(),
            self.required_coverage_threshold, ', '.join(runs)
        )
        if non_useable_runs:
            s += '(non useable run elements in %s)' % ', '.join(non_useable_runs)
        else:
            s += '(no non useable run elements)'
        return s

    def start(self):
        super().start()
        self.lims_ntf.start_sample_pipeline()

    def succeed(self):
        super().succeed()
        self.lims_ntf.assign_next_and_advance_step()

    def fail(self, exit_status):
        super().fail(exit_status)
        self.lims_ntf.remove_sample_from_workflow()

    @property
    def pipeline(self):
        analysis_type = self.lims_sample_info.get('Analysis Type')

        if self.species is None:
            raise AnalysisDriverError('No species information found in the LIMS for ' + self.name)

        elif self.species == 'Homo sapiens':
            if 'Variant Calling gatk4' in analysis_type:
                cls = human_variant_calling_gatk4.HumanVarCallingGATK4
            else:
                cls = bcbio.BCBioVarCalling
        elif analysis_type in ['Variant Calling gatk4']:
            cls = variant_calling_gatk4.VarCallingGATK4
        elif analysis_type in ['Variant Calling', 'Variant Calling gatk', 'Variant Calling gatk3']:
            cls = variant_calling.VarCalling
        elif analysis_type in ['QC GATK3']:
            # This is unlikely to be used in production but allows us to trigger the GATK3 QC pipeline when needed
            cls = qc.QC
        else:
            cls = qc_gatk4.QCGATK4

        return cls(self)


class ProjectDataset(Dataset):
    type = 'project'
    endpoint = 'projects'
    id_field = 'project_id'

    def __init__(self, name, most_recent_proc=None):
        super().__init__(name, most_recent_proc)
        self._number_of_samples = None
        self._samples_processed = None
        self._lims_project_info = None
        self._lims_samples_info = None
        self._sample_datasets = {}
        self._species = None
        self._genome_version = None

    def _is_ready(self):
        return 0 < self.number_of_samples <= len(self.samples_processed)

    @property
    def data_source(self):
        return [s['sample_id'] for s in self.samples_processed]

    @property
    def samples_processed(self):
        if not self._samples_processed:
            self._samples_processed = rest_communication.get_documents(
                'samples',
                where={'project_id': self.name, 'aggregated.most_recent_proc.status': 'finished'},
                all_pages=True
            )
        return self._samples_processed

    def sample_dataset(self, sample_id):
        if sample_id not in self._sample_datasets:
            self._sample_datasets[sample_id] = SampleDataset(sample_id)
        return self._sample_datasets[sample_id]

    @property
    def lims_project_info(self):
        if self._lims_project_info is None:
            self._lims_project_info = rest_communication.get_document(
                'lims/project_info', match={'project_id': self.name}
            )
        return self._lims_project_info

    def get_processed_gvcfs(self):
        project_source = os.path.join(cfg.query('project', 'input_dir'), self.name)
        gvcf_generating_pipelines = ('bcbio', 'variant_calling_gatk4', 'human_variant_calling_gatk4')
        gvcf_files = []
        for sample in self.samples_processed:
            if query_dict(sample, 'aggregated.most_recent_proc.pipeline_used.name') in gvcf_generating_pipelines:
                gvcf_file = find_file(project_source, sample['sample_id'], sample['user_sample_id'] + '.g.vcf.gz')
                if not gvcf_file:
                    raise AnalysisDriverError(
                        'Unable to find gVCF file for sample %s in %s' % (sample['sample_id'], project_source)
                    )
                gvcf_files.append(gvcf_file)

        return gvcf_files

    @property
    def number_of_samples(self):
        if not self._number_of_samples:
            self._number_of_samples = int(self.lims_project_info.get('nb_quoted_samples', -1))
        return self._number_of_samples

    @property
    def species(self):
        if self._species is None:
            s = set()
            for sample in self.samples_processed:
                s.add(sample.get('species_name') or self.sample_dataset(sample['sample_id']).species)
            if len(s) != 1:
                raise AnalysisDriverError('Unexpected number of species (%s) in this project' % ', '.join(s))
            self._species = s.pop()
        return self._species

    @property
    def genome_version(self):
        if self._genome_version is None:
            self._genome_version = self.sample_dataset(self.samples_processed[0]['sample_id']).genome_version
        return self._genome_version

    @property
    def reference_genome(self):
        return self.sample_dataset(self.samples_processed[0]['sample_id']).reference_genome

    def __str__(self):
        return '%s  (%s samples / %s) ' % (super().__str__(), len(self.samples_processed), self.number_of_samples)

    @property
    def pipeline(self):
        return projects.Project(self)


class MostRecentProc:
    def __init__(self, dataset, initial_content=None):
        self.lock = Lock()
        self.dataset = dataset
        self.proc_id = None
        if initial_content:
            self._entity = initial_content.copy()
        else:
            self._entity = None

    def retrieve_entity(self):
        procs = rest_communication.get_documents(
            'analysis_driver_procs',
            where={'dataset_type': self.dataset.type, 'dataset_name': self.dataset.name},
            sort='-_created'
        )
        if procs:
            # remove the private (starting with '_') fields from the dict
            self._entity = {k: v for k, v in procs[0].items() if not k.startswith('_')}
            self.proc_id = self.entity['proc_id']
        else:
            self._entity = {}

    @property
    def entity(self):
        if self._entity is None:
            self.retrieve_entity()
        return self._entity

    def initialise_entity(self):
        self.proc_id = '%s_%s_%s' % (self.dataset.type, self.dataset.name, now())
        entity = {
            'proc_id': self.proc_id,
            'dataset_type': self.dataset.type,
            'dataset_name': self.dataset.name
        }

        if self.dataset.data_source:
            entity.update({'data_source': self.dataset.data_source})

        rest_communication.post_entry('analysis_driver_procs', entity)
        rest_communication.patch_entry(
            self.dataset.endpoint,
            {'analysis_driver_procs': [self.proc_id]},
            id_field=self.dataset.id_field,
            element_id=self.dataset.name,
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
        if not self.entity:  # initialise self._entity with database content, or {} if none available
            self.initialise_entity()  # if self._entity == {}, then initialise and push to database

        self.entity.update(kwargs)
        self.sync()

    def change_status(self, status):
        with self.lock:
            self.update_entity(status=status)

    def start(self):
        with self.lock:
            payload = {
                'pipeline_used': {
                    'name': self.dataset.pipeline.name,
                    'toolset_type': toolset.type,
                    'toolset_version': toolset.version
                }
            }
            if self.get('status') != DATASET_RESUME:
                payload['start_date'] = now()

            if self.dataset.type == 'sample':
                payload['genome_used'] = self.dataset.genome_version

            self.update_entity(status=DATASET_PROCESSING, pid=os.getpid())
            rest_communication.patch_entry('analysis_driver_procs', payload, 'proc_id', self.proc_id)

    def finish(self, status):
        with self.lock:
            self.update_entity(status=status, pid=None, end_date=now())

    def start_stage(self, stage_name):
        with self.lock:
            doc = rest_communication.get_document(
                'analysis_driver_stages',
                where={'stage_id': self._stage_id(stage_name)}
            )
            if doc:
                rest_communication.patch_entry(
                    'analysis_driver_stages',
                    {'date_started': now(), 'date_finished': None, 'exit_status': None},
                    'stage_id', self._stage_id(stage_name)
                )
            else:
                rest_communication.post_entry(
                    'analysis_driver_stages',
                    {'stage_id': self._stage_id(stage_name), 'date_started': now(),
                     'stage_name': stage_name, 'analysis_driver_proc': self.proc_id}
                )
                self.retrieve_entity()
                stages = self.entity.get('stages', [])
                stages.append(self._stage_id(stage_name))
                self.update_entity(stages=stages)

    def end_stage(self, stage_name, exit_status=0):
        rest_communication.patch_entry(
            'analysis_driver_stages',
            {'date_finished': now(), 'exit_status': exit_status}, 'stage_id', self._stage_id(stage_name)
        )

    def get(self, key, ret_default=None):
        return self.entity.get(key, ret_default)

    def _stage_id(self, stage_name):
        return self.proc_id + '_' + stage_name
