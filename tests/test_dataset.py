import os
import time
import pytest
from sys import modules
from unittest.mock import patch, Mock, PropertyMock
from egcg_core import constants as c
from integration_tests.mocked_data import MockedSample
from tests.test_analysisdriver import TestAnalysisDriver, NamedMock
from analysis_driver.exceptions import AnalysisDriverError, SequencingRunError
from analysis_driver.dataset import Dataset, RunDataset, SampleDataset, ProjectDataset, MostRecentProc
from analysis_driver.tool_versioning import Toolset

ppath = 'analysis_driver.dataset.'
fake_proc = {'proc_id': 'a_proc_id', 'dataset_type': 'test', 'dataset_name': 'test'}

patched_patch = patch(ppath + 'rest_communication.patch_entry')
patched_post = patch(ppath + 'rest_communication.post_entry')
patched_update = patch(ppath + 'MostRecentProc.update_entity')


def patched_get_doc(content=None):
    return patch(ppath + 'rest_communication.get_document', return_value=content or fake_proc)


def patched_get_docs(content=None):
    return patch(ppath + 'rest_communication.get_documents', return_value=content or [fake_proc])


def patched_datetime(time='now'):
    return patch(ppath + 'now', return_value=time)


class FakeProc(Mock):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.entity = {}

    def get(self, key, ret_default=None):
        return self.entity.get(key, ret_default)


class TestDataset(TestAnalysisDriver):
    base_dir = os.path.join(TestAnalysisDriver.assets_path, 'dataset_scanner')

    def setUp(self):
        self.dataset = Dataset('a_dataset')
        self.dataset._ntf = Mock()
        self.dataset.most_recent_proc = FakeProc()

    def test_dataset_status(self):
        for s in (c.DATASET_REPROCESS, None):
            self.dataset.most_recent_proc.entity['status'] = s

            with patch.object(self.dataset.__class__, '_is_ready') as mocked_ready:
                mocked_ready.return_value = True
                assert self.dataset.dataset_status == c.DATASET_READY
                mocked_ready.return_value = False
                assert self.dataset.dataset_status == c.DATASET_NEW

        self.dataset.most_recent_proc.entity['status'] = 'a_status'
        assert self.dataset.dataset_status == 'a_status'

    @patched_get_docs([{'stage_name': 'this'}, {'stage_name': 'that'}])
    def test_running_stages(self, mocked_get):
        self.dataset.most_recent_proc.entity['proc_id'] = 'a_proc_id'
        assert self.dataset.running_stages == ['this', 'that']
        mocked_get.assert_called_with(
            'analysis_driver_stages',
            all_pages=True,
            where={'analysis_driver_proc': 'a_proc_id', 'date_finished': None}
        )

    def test_start(self):
        self.dataset.most_recent_proc.entity['status'] = c.DATASET_PROCESSING
        with self.assertRaises(AssertionError):
            self.dataset.start()

        self.dataset.most_recent_proc.retrieve_entity.reset_mock()
        self.dataset.most_recent_proc.entity['status'] = c.DATASET_RESUME
        self.dataset.start()

        # 1 from start, 1 from start -> _assert_status
        assert self.dataset.most_recent_proc.retrieve_entity.call_count == 2
        assert self.dataset.most_recent_proc.start.call_count == 1
        assert self.dataset.ntf.start_pipeline.call_count == 1

        with patch.object(self.dataset.__class__, 'initialise_entity') as mocked_init:
            self.dataset.most_recent_proc.entity['status'] = c.DATASET_READY
            self.dataset.start()
            assert mocked_init.call_count == 1

    def test_succeed(self):
        self.dataset.most_recent_proc.entity['status'] = c.DATASET_PROCESSING
        self.dataset.succeed()
        self.dataset.most_recent_proc.finish.assert_called_with(c.DATASET_PROCESSED_SUCCESS)
        self.dataset.ntf.end_pipeline.assert_called_with(0)

    def test_fail(self):
        self.dataset.most_recent_proc.entity['status'] = c.DATASET_PROCESSING
        self.dataset.fail(9)
        self.dataset.most_recent_proc.finish.assert_called_with(c.DATASET_PROCESSED_FAIL)
        self.dataset.ntf.end_pipeline.assert_called_with(9)

    def test_abort(self):
        self.dataset.abort()
        self.dataset.most_recent_proc.finish.assert_called_with(c.DATASET_ABORTED)

    def test_resume(self):
        with patch.object(self.dataset.__class__, 'terminate') as mocked_term:
            self.dataset.resume()
            assert mocked_term.call_count == 1
            self.dataset.most_recent_proc.change_status.assert_called_with(c.DATASET_RESUME)

    def test_reset(self):
        with patch.object(self.dataset.__class__, 'terminate') as mocked_term:
            self.dataset.reset()
            assert mocked_term.call_count == 1
            self.dataset.most_recent_proc.change_status.assert_called_with(c.DATASET_REPROCESS)

    def test_skip(self):
        self.dataset.skip()
        self.dataset.most_recent_proc.finish.assert_called_with(c.DATASET_PROCESSED_SUCCESS)

    def test_force(self):
        self.dataset.force()
        self.dataset.most_recent_proc.change_status.assert_called_with(c.DATASET_FORCE_READY)

    @patch('os.kill')
    def test_terminate(self, mocked_kill):
        self.dataset.most_recent_proc.entity = {'pid': 1337}

        self.dataset.__dict__['_logger'] = Mock()
        with patch.object(self.dataset.__class__, '_pid_valid', return_value=False) as mocked_valid:
            self.dataset._terminate(1)
            mocked_valid.assert_called_with(1337)
            self.dataset._logger.debug.assert_called_with('Attempted to terminate invalid pid %s with signal %s', 1337, 1)

            with patch.object(self.dataset.__class__, '_pid_running', return_value=False):
                mocked_valid.return_value = True
                self.dataset._terminate(2)
                mocked_kill.assert_called_with(1337, 2)
                self.dataset._logger.info.assert_called_with(
                    'Terminated pid %s for %s %s with signal %s',
                    1337, self.dataset.type, 'a_dataset', 2
                )

    def test_pid_valid(self):
        cmdlinefile = os.path.join(self.assets_path, 'example.pid')
        if os.path.isfile(cmdlinefile):
            os.remove(cmdlinefile)

        with patch('os.path.join', return_value=cmdlinefile):
            assert self.dataset._pid_valid(1337) is False

            with open(cmdlinefile, 'w') as f:
                f.write(modules['__main__'].__file__ + '\n')

            assert self.dataset._pid_valid(1337) is True

    def test_start_stage(self):
        self.dataset.start_stage('a_stage')
        self.dataset.most_recent_proc.start_stage.assert_called_with('a_stage')
        self.dataset.ntf.start_stage.assert_called_with('a_stage')

    def test_end_stage(self):
        self.dataset.end_stage('a_stage')
        self.dataset.most_recent_proc.end_stage.assert_called_with('a_stage', 0)
        self.dataset.ntf.end_stage.assert_called_with('a_stage', 0)

    def test_assert_status(self):
        with patch.object(self.dataset.__class__, 'dataset_status', new=c.DATASET_NEW):
            self.dataset._assert_status(c.DATASET_NEW)
            with self.assertRaises(AssertionError) as e:
                self.dataset._assert_status(c.DATASET_READY)
            assert str(e.exception) == 'Status assertion failed on a new %s' % self.dataset.type
            assert self.dataset.most_recent_proc.retrieve_entity.call_count == 2

    def test_report(self):
        self.dataset.most_recent_proc.entity = {'pid': 1337}
        with patch.object(self.dataset.__class__, 'running_stages', new=['this', 'that']):
            assert self.dataset.report() == 'a_dataset (1337) -- this, that'

    def test_raise_exceptions(self):
        self.dataset.register_exception(Mock(stage_name='task1'), SequencingRunError('RunErrored'))
        # Give time for the exception to be pickled.
        # From the documentation:
        # 'After putting an object on an empty queue there may be an infinitesimal delay before the queue’s empty()
        # method returns False and get_nowait() can return without raising queue.Empty.'
        time.sleep(.1)
        with self.assertRaises(SequencingRunError) as e:
            self.dataset.raise_exceptions()

        assert str(e.exception) == 'RunErrored'

    def test_resolve_pipeline_and_toolset(self):
        self.dataset._pipeline_instruction = Mock(
            return_value={'toolset_type': 'a_toolset', 'toolset_version': 3, 'name': 'qc'}
        )
        with patch(ppath + 'toolset') as mocked_toolset:
            self.dataset.resolve_pipeline_and_toolset()

        mocked_toolset.select_type.assert_called_with('a_toolset')
        mocked_toolset.select_version.assert_called_with(3)
        assert self.dataset.pipeline.name == 'qc'

    @patch(ppath + 'toolset', new=Toolset())
    def test_pipeline_instruction(self):
        self.dataset._default_pipeline = Mock(return_value='qc')
        assert self.dataset._pipeline_instruction() == {
            'name': 'qc',
            'toolset_type': 'non_human_sample_processing',
            'toolset_version': 0
        }

mocked_lane_user_prep_artifact1 = NamedMock(real_name='art1', reagent_labels=[], samples=[MockedSample(real_name='sample1')])
mocked_lane_user_prep_artifact2 = NamedMock(real_name='art2', reagent_labels=[], samples=[MockedSample(real_name='sample2')])

mocked_lane_artifact1 = NamedMock(real_name='art1', reagent_labels=['D701-D502 (ATTACTCG-ATAGAGGC)'], samples=[MockedSample(real_name='sample1')])
mocked_lane_artifact2 = NamedMock(real_name='art2', reagent_labels=['D702-D502 (TCCGGAGA-ATAGAGGC)'], samples=[MockedSample(real_name='sample2')])
mocked_lane_artifact3 = NamedMock(real_name='art3', reagent_labels=['D703-D502 (CGCTCATT-ATAGAGGC)'], samples=[MockedSample(real_name='sample3')])
mocked_lane_artifact4 = NamedMock(real_name='art4', reagent_labels=['D704-D502 (GAGATTCC-ATAGAGGC)'], samples=[MockedSample(real_name='sample4')])
mocked_lane_artifact5 = NamedMock(real_name='art5', reagent_labels=['D705-D502 (ATTCAGAA-ATAGAGGC)'], samples=[MockedSample(real_name='sample5')])
mocked_lane_artifact6 = NamedMock(real_name='art6', reagent_labels=['D706-D502 (GAATTCGT-ATAGAGGC)'], samples=[MockedSample(real_name='sample6')])
mocked_lane_artifact_pool = NamedMock(real_name='artpool', reagent_labels=[
    'D703-D502 (CGCTCATT-ATAGAGGC)',
    'D704-D502 (GAGATTCC-ATAGAGGC)',
    'D705-D502 (ATTCAGAA-ATAGAGGC)',
    'D706-D502 (GAATTCGT-ATAGAGGC)',
], samples=[
    'sample3',
    'sample4',
    'sample5',
    'sample6'
], input_artifact_list=Mock(return_value=[
    mocked_lane_artifact3,
    mocked_lane_artifact4,
    mocked_lane_artifact5,
    mocked_lane_artifact6
]), parent_process=Mock(type=NamedMock(real_name='Create PDP Pool')))

mocked_flowcell_non_pooling = Mock(placements={
    '1:1': mocked_lane_artifact1,
    '2:1': mocked_lane_artifact2,
    '3:1': mocked_lane_artifact1,
    '4:1': mocked_lane_artifact2,
    '5:1': mocked_lane_artifact1,
    '6:1': mocked_lane_artifact2,
    '7:1': mocked_lane_artifact1,
    '8:1': mocked_lane_artifact2
})

mocked_flowcell_pooling = Mock(placements={
    '1:1': mocked_lane_artifact_pool,
    '2:1': mocked_lane_artifact_pool,
    '3:1': mocked_lane_artifact_pool,
    '4:1': mocked_lane_artifact_pool,
    '5:1': mocked_lane_artifact_pool,
    '6:1': mocked_lane_artifact_pool,
    '7:1': mocked_lane_artifact_pool,
    '8:1': mocked_lane_artifact_pool
})

mocked_flowcell_user_prepared = Mock(placements={
    '1:1': mocked_lane_user_prep_artifact1,
    '2:1': mocked_lane_user_prep_artifact2,
    '3:1': mocked_lane_user_prep_artifact1,
    '4:1': mocked_lane_user_prep_artifact2,
    '5:1': mocked_lane_user_prep_artifact1,
    '6:1': mocked_lane_user_prep_artifact2,
    '7:1': mocked_lane_user_prep_artifact1,
    '8:1': mocked_lane_user_prep_artifact2
})


class TestRunDataset(TestDataset):
    def setUp(self):
        self.dataset = RunDataset('a_dataset')
        self.dataset._ntf = Mock()
        self.dataset.most_recent_proc = FakeProc()

    @patched_post
    @patched_get_doc('a_run')
    def test_initialise_entity(self, mocked_get, mocked_post):
        self.dataset.initialise_entity()
        assert mocked_post.call_count == 0

        mocked_get.return_value = None
        self.dataset.initialise_entity()
        assert mocked_post.call_count == 1

    def test_is_ready(self):
        assert self.dataset._is_ready()

    def test_find_pooling_step(self):
        input_artifacts = ['an_input_artifact', 'another_input_artifact']

        def mocked_artifact(input_artifacts):
            return Mock(
                input_artifact_list=Mock(return_value=input_artifacts),
                parent_process=Mock(type=NamedMock(real_name='Create PDP Pool'))
            )

        art3 = mocked_artifact(input_artifacts)
        art2 = mocked_artifact([art3])
        art1 = mocked_artifact([art2])

        assert self.dataset._find_pooling_step_for_artifact(art1) == input_artifacts
        with self.assertRaises(ValueError) as e:
            self.dataset._find_pooling_step_for_artifact(art1, 'Create PDP Pool', 2)
        assert str(e.exception) == 'Cannot find pooling step after 2 iterations'

        art3.parent_process.type.real_name = 'Not Create PDP Pool'
        with self.assertRaises(ValueError) as e:
            self.dataset._find_pooling_step_for_artifact(art1, 'Create PDP Pool')
        assert str(e.exception) == 'Mismatching step name: Create PDP Pool != Not Create PDP Pool'

    def test_run_elements_from_lims(self):

        def patched_run(fake_flowcell):
            return patch(
                'analysis_driver.dataset.clarity.get_run',
                return_value=Mock(
                    parent_processes=Mock(return_value=[Mock(output_containers=Mock(return_value=[fake_flowcell]))])
                )
            )

        def patched_barcodes(has_barcodes):
            return patch.object(RunDataset, 'has_barcodes', new_callable=PropertyMock(return_value=has_barcodes))

        d = RunDataset('test_dataset')
        with patched_run(mocked_flowcell_non_pooling), patched_barcodes(False):
            run_elements = d._run_elements_from_lims()
            assert len(set(r[c.ELEMENT_PROJECT_ID] for r in run_elements)) == 1
            assert len(set(r[c.ELEMENT_SAMPLE_INTERNAL_ID] for r in run_elements)) == 2
            barcodes_len = set(len(r[c.ELEMENT_BARCODE]) for r in run_elements)
            assert len(barcodes_len) == 1
            assert barcodes_len.pop() == 0

        d = RunDataset('test_dataset')
        with patched_run(mocked_flowcell_pooling), patched_barcodes(True):
            run_elements = d._run_elements_from_lims()
            assert len(set(r[c.ELEMENT_PROJECT_ID] for r in run_elements)) == 1
            assert len(set(r[c.ELEMENT_SAMPLE_INTERNAL_ID] for r in run_elements)) == 4
            barcodes_len = set(len(r[c.ELEMENT_BARCODE]) for r in run_elements)
            assert len(barcodes_len) == 1
            assert barcodes_len.pop() == 8

        d = RunDataset('test_dataset')
        with patched_run(mocked_flowcell_user_prepared), patched_barcodes(False):
            run_elements = d._run_elements_from_lims()
            assert len(set(r[c.ELEMENT_PROJECT_ID] for r in run_elements)) == 1
            assert len(set(r[c.ELEMENT_SAMPLE_INTERNAL_ID] for r in run_elements)) == 2
            barcodes_len = set(len(r[c.ELEMENT_BARCODE]) for r in run_elements)
            assert len(barcodes_len) == 1
            assert barcodes_len.pop() == 0

    @patch('builtins.open')
    def test_generate_samplesheet(self, mocked_open):
        fake_run_elements = [
            {'lane': '1', 'sample_id': 'sample_1', 'library_id': 'lib_1', 'project_id': 'p', 'barcode': 'ATGC'},
            {'lane': '2', 'sample_id': 'sample_2', 'library_id': 'lib_2', 'project_id': 'p', 'barcode': 'CTGA'}
        ]
        exp = [
            '[Header]', 'Date, now', 'Workflow, Generate FASTQ Only', '',
            '[Settings]', 'Adapter, AGATCGGAAGAGCACACGTCTGAACTCCAGTCA',
            'AdapterRead2, AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT', '', '[Data]',
            'Lane,Sample_ID,Sample_Name,Sample_Project,index',
            '1,sample_1,lib_1,p,ATGC', '2,sample_2,lib_2,p,CTGA', ''
        ]

        self.dataset._run_elements = fake_run_elements
        with patched_datetime():
            self.dataset._generate_samplesheet('a_samplesheet')
            mocked_open.return_value.__enter__.return_value.write.assert_called_with('\n'.join(exp))

    @patch('analysis_driver.dataset.clarity.get_run')
    def test_is_sequencing(self, mocked_get_run):
        mocked_get_run.return_value = None
        with patch.object(RunDataset, 'warning') as mocked_log:
            assert self.dataset.is_sequencing()
            mocked_log.assert_called_with('Run %s not found in the LIMS', 'a_dataset')
        mocked_get_run.return_value = Mock(udf={'Run Status': 'RunStarted'})
        assert self.dataset.is_sequencing()
        mocked_get_run.return_value.get.assert_called_with(force=True)


class TestSampleDataset(TestDataset):
    def setUp(self):
        self.dataset = SampleDataset('a_dataset')
        self.dataset._ntf = Mock()
        self.dataset._lims_ntf = Mock()
        self.dataset.most_recent_proc = FakeProc()

    def test_amount_data(self):
        self.dataset._sample = {'aggregated': {'clean_yield_in_gb': 15}}
        assert self.dataset._amount_data() == 15000000000

    def test_data_threshold(self):
        self.dataset._sample = {}
        with self.assertRaises(AnalysisDriverError) as e:
            _ = self.dataset.data_threshold
        assert str(e.exception) == 'Could not find data threshold for a_dataset'

        self.dataset._sample = {'required_yield': 1000000000}
        assert self.dataset.data_threshold == 1000000000

    def test_no_data_threshold(self):
        self.dataset._sample = {}
        with pytest.raises(AnalysisDriverError) as e:
            _ = self.dataset.data_threshold
        assert 'Could not find data threshold' in str(e)

    @patch(ppath + 'toolset', new=Toolset())
    @patch(ppath + 'SampleDataset.project_id', new='a_project')
    def test_pipeline_instruction(self):
        self.dataset._default_pipeline = Mock(return_value='qc')

        with patched_get_doc({'sample_pipeline': 'some_data'}):
            assert self.dataset._pipeline_instruction() == 'some_data'

        exp = {
            'name': 'qc',
            'toolset_type': 'non_human_sample_processing',
            'toolset_version': 0
        }

        with patched_get_doc(None), patched_patch as mocked_patch:
            assert self.dataset._pipeline_instruction() == exp
            mocked_patch.assert_called_with('projects', {'sample_pipeline': exp}, 'project_id', 'a_project')

    def test_is_ready(self):
        self.dataset._sample = {'aggregated': {'clean_pc_q30': 85, 'clean_yield_in_gb': 1.5}}
        self.dataset._data_threshold = 20000000000
        assert not self.dataset._is_ready()
        self.dataset._data_threshold = 100
        assert self.dataset._is_ready()

    def test_report(self):
        self.dataset.most_recent_proc.entity = {'pid': 1337}
        self.dataset._sample = {'aggregated': {'run_ids': ['other', 'another'], 'clean_yield_in_gb': 1}}
        self.dataset._non_useable_run_elements = []
        self.dataset._data_threshold = 2000000000
        exp = 'a_dataset (1337) -- this, that  (1000000000 / 2000000000  from other, another)'

        with patch.object(SampleDataset, 'running_stages', new=['this', 'that']):
            assert self.dataset.report() == exp + ' (no non useable run elements)'
            self.dataset._non_useable_run_elements = [{'run_id': 'more'}]
            assert self.dataset.report() == exp + ' (non useable run elements in more)'

    def test_start(self):
        super().test_start()
        assert self.dataset.lims_ntf.start_sample_pipeline.call_count == 2

    def test_succeed(self):
        super().test_succeed()
        assert self.dataset.lims_ntf.assign_next_and_advance_step.call_count == 1

    def test_fail(self):
        super().test_fail()
        assert self.dataset.lims_ntf.remove_sample_from_workflow.call_count == 1

    def test_data_source(self):
        self.dataset._run_elements = [{'run_element_id': 'run_element_1'}, {'run_element_id': 'run_element_2'}]
        assert self.dataset.data_source == ['run_element_1', 'run_element_2']

    def test_default_pipeline(self):
        self.dataset._species = 'Homo sapiens'
        assert self.dataset._default_pipeline() == 'bcbio'

        self.dataset._species = 'Thingius thingy'
        patched_get_sample = patch(
            'analysis_driver.dataset.clarity.get_sample',
            return_value=Mock(udf={'Analysis Type': 'Variant Calling'})
        )
        with patched_get_sample as mocked_sample:
            assert self.dataset._default_pipeline() == 'variant_calling'
            mocked_sample.return_value.udf = {}
            assert self.dataset._default_pipeline() == 'qc'


class TestProjectDataset(TestDataset):
    def setUp(self):
        self.dataset = ProjectDataset('a_dataset')
        self.dataset._ntf = Mock()
        self.dataset.most_recent_proc = FakeProc()

    def test_is_ready(self):
        self.dataset._number_of_samples = 3
        self.dataset._samples_processed = ['this', 'that']
        assert not self.dataset._is_ready()
        self.dataset._samples_processed.append('other')
        assert self.dataset._is_ready()

    @patch('analysis_driver.dataset.clarity.get_project', return_value=Mock(udf={'Number of Quoted Samples': 3}))
    def test_number_of_samples(self, mocked_get_project):
        assert self.dataset.number_of_samples == 3
        self.dataset._number_of_samples = None
        mocked_get_project.return_value.udf['Number of Quoted Samples'] = 0
        assert self.dataset.number_of_samples == -1
        self.dataset._number_of_samples = None
        mocked_get_project.return_value = None
        assert self.dataset.number_of_samples == -1

    @patch('analysis_driver.dataset.clarity.get_species_from_sample', return_value='Thingius thingy')
    def test_species(self, mocked_get_species):
        self.dataset._samples_processed = [
            {'sample_id': 'a_sample', 'species_name': 'Thingius thingy'},
            {'sample_id': 'another_sample'}
        ]
        assert self.dataset.species == 'Thingius thingy'

    @patch('analysis_driver.dataset.clarity.get_sample', return_value=Mock(udf={'Genome Version': 'a_genome'}))
    def test_genome_version(self, mocked_get_sample):
        self.dataset._samples_processed = [{'sample_id': 'a_sample'}]
        self.dataset._species = 'Homo sapiens'
        assert self.dataset.genome_version == 'a_genome'

        self.dataset._genome_version = None
        mocked_get_sample.return_value.udf = {}
        assert self.dataset.genome_version == 'genome_version'  # from config


class TestMostRecentProc(TestAnalysisDriver):
    def setUp(self):
        with patched_get_docs():
            self.proc = MostRecentProc(
                NamedMock(
                    type='test',
                    real_name='test',
                    pipeline_type='some kind of variant calling',
                    data_source=None,
                    endpoint='tests',
                    id_field='test_id'
                )
            )
        self.proc._entity = fake_proc.copy()
        self.proc.proc_id = 'a_proc_id'

    def test_rest_entity_pre_existing(self):
        assert self.proc.entity == fake_proc

    @patched_get_docs()
    @patch(ppath + 'MostRecentProc.initialise_entity')
    def test_rest_entity_not_pre_existing(self, mocked_initialise, mocked_get):
        self.proc._entity = None
        with patched_datetime():
            x = self.proc.entity
            assert x == fake_proc
            mocked_get.assert_called_with(
                'analysis_driver_procs',
                where={'dataset_type': 'test', 'dataset_name': 'test'},
                sort='-_created'
            )

        self.proc._entity = None
        mocked_get.return_value = []
        y = self.proc.entity
        assert y == {}
        mocked_get.assert_called_with(
            'analysis_driver_procs',
            where={'dataset_type': 'test', 'dataset_name': 'test'},
            sort='-_created'
        )
        assert mocked_initialise.call_count == 0

    @patched_post
    @patched_get_docs()
    @patched_patch
    def test_data_source(self, mocked_patch, mocked_get, mocked_post):
        sample_proc = MostRecentProc(
            NamedMock(
                type='sample',
                real_name='sample1',
                pipeline_type='some kind of variant calling',
                data_source=['run_element1', 'run_element2']
            )
        )
        sample_proc._entity = None
        _ = sample_proc.entity
        with patched_datetime():
            sample_proc.initialise_entity()
            mocked_post.assert_called_with(
                'analysis_driver_procs',
                {
                    'dataset_type': 'sample',
                    'dataset_name': 'sample1',
                    'data_source': ['run_element1', 'run_element2'],
                    'proc_id': 'sample_sample1_now'
                }
            )

        project_proc = MostRecentProc(
            NamedMock(
                type='project',
                real_name='project1',
                pipeline_type='project pipeline',
                data_source=['sample_1', 'sample_2']
            )
        )
        project_proc._entity = None
        _ = project_proc.entity
        with patched_datetime():
            project_proc.initialise_entity()
            mocked_post.assert_called_with(
                'analysis_driver_procs',
                {
                    'data_source': ['sample_1', 'sample_2'],
                    'dataset_name': 'project1',
                    'proc_id': 'project_project1_now',
                    'dataset_type': 'project'
                }
            )

    @patched_post
    @patched_patch
    def test_initialise_entity(self, mocked_patch, mocked_post):
        self.proc._entity = None
        with patched_datetime():
            self.proc.initialise_entity()
        mocked_post.assert_called_with(
            'analysis_driver_procs',
            {'proc_id': 'test_test_now', 'dataset_type': 'test', 'dataset_name': 'test'}
        )
        mocked_patch.assert_called_with(
            'tests',
            {'analysis_driver_procs': ['test_test_now']},
            id_field='test_id',
            element_id='test',
            update_lists=['analysis_driver_procs']
        )

        with patched_datetime('then'):
            self.proc.initialise_entity()
        mocked_post.assert_called_with(
            'analysis_driver_procs',
            {'proc_id': 'test_test_then', 'dataset_type': 'test', 'dataset_name': 'test'}
        )
        mocked_patch.assert_called_with(
            'tests',
            {'analysis_driver_procs': ['test_test_then']},
            id_field='test_id',
            element_id='test',
            update_lists=['analysis_driver_procs']
        )

    @patched_patch
    def test_sync(self, mocked_patch):
        self.proc.sync()
        assert self.proc._entity == {
            'proc_id': 'a_proc_id',
            'dataset_type': 'test',
            'dataset_name': 'test'
        }

        self.proc.entity.update({'this': 'that'})
        self.proc.sync()
        mocked_patch.assert_called_with(
            'analysis_driver_procs',
            {'this': 'that', 'dataset_type': 'test', 'dataset_name': 'test'},
            id_field='proc_id',
            element_id='a_proc_id'
        )
        assert self.proc._entity == {
            'proc_id': 'a_proc_id',
            'dataset_type': 'test',
            'dataset_name': 'test',
            'this': 'that'
        }

    @patch(ppath + 'MostRecentProc.sync')
    def test_update_entity(self, mocked_sync):
        self.proc.update_entity(other='another')
        assert self.proc.entity == {'proc_id': 'a_proc_id', 'dataset_type': 'test',
                                    'dataset_name': 'test', 'other': 'another'}
        mocked_sync.assert_called()

    @patch(ppath + 'now', return_value='now')
    @patch(ppath + 'toolset', new=Mock(type='some kind of toolset', version=3))
    @patched_patch
    @patched_update
    def test_start(self, mocked_update, mocked_patch, mocked_now):
        self.proc.dataset.pipeline = NamedMock(real_name='some kind of variant calling')
        with patch('os.getpid', return_value=1):
            self.proc.start()
        mocked_update.assert_called_with(status=c.DATASET_PROCESSING, pid=1)
        mocked_patch.assert_called_with(
            'analysis_driver_procs',
            {
                'start_date': 'now',
                'pipeline_used': {
                    'name': 'some kind of variant calling',
                    'toolset_type': 'some kind of toolset',
                    'toolset_version': 3
                }
            },
            'proc_id',
            'a_proc_id'
        )

    @patched_update
    def test_finish(self, mocked_update):
        with patched_datetime():
            self.proc.finish(c.DATASET_PROCESSED_SUCCESS)
        mocked_update.assert_called_with(status=c.DATASET_PROCESSED_SUCCESS, pid=None, end_date='now')

    @patched_update
    @patched_post
    @patched_patch
    @patch(ppath + 'MostRecentProc.retrieve_entity')
    @patch('analysis_driver.dataset.rest_communication.get_document')
    def test_start_stage(self, mocked_get_doc, mocked_retrieve, mocked_patch, mocked_post, mocked_update):
        with patched_datetime('then'):
            mocked_get_doc.return_value = True
            self.proc.start_stage('stage_1')
            assert mocked_update.call_count == 0
            mocked_patch.assert_called_with(
                'analysis_driver_stages',
                {'date_started': 'then', 'date_finished': None, 'exit_status': None},
                'stage_id', 'a_proc_id_stage_1'
            )

            mocked_get_doc.return_value = None
            self.proc.start_stage('stage_1')
            mocked_update.assert_called_with(stages=['a_proc_id_stage_1'])
            mocked_post.assert_called_with(
                'analysis_driver_stages',
                {'stage_id': 'a_proc_id_stage_1', 'date_started': 'then', 'stage_name': 'stage_1', 'analysis_driver_proc': 'a_proc_id'}
            )

        self.proc.entity['stages'] = ['a_proc_id_stage_1']
        with patched_datetime('later'):
            self.proc.start_stage('stage_2')

        mocked_update.assert_called_with(stages=['a_proc_id_stage_1', 'a_proc_id_stage_2'])
        mocked_post.assert_called_with(
            'analysis_driver_stages',
            {'stage_id': 'a_proc_id_stage_2', 'date_started': 'later', 'stage_name': 'stage_2', 'analysis_driver_proc': 'a_proc_id'}
        )

    @patched_patch
    @patched_datetime()
    def test_end_stage(self, mocked_now, mocked_patch):
        self.proc.end_stage('stage_1')
        mocked_patch.assert_called_with(
            'analysis_driver_stages', {'date_finished': 'now', 'exit_status': 0}, 'stage_id', 'a_proc_id_stage_1'
        )

    def test_stage_id(self):
        self.proc.proc_id = 'test_proc'
        assert self.proc._stage_id('stage_name') == 'test_proc_stage_name'
