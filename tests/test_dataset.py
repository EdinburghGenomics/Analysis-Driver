import os
import signal
import time

import pytest
from sys import modules
from unittest.mock import patch, Mock, PropertyMock, call
from egcg_core import constants as c
from integration_tests.mocked_data import MockedSample, MockedRunProcess
from tests.test_analysisdriver import TestAnalysisDriver, NamedMock
from analysis_driver.exceptions import AnalysisDriverError, SequencingRunError
from analysis_driver.dataset import Dataset, RunDataset, SampleDataset, ProjectDataset, MostRecentProc
from analysis_driver.tool_versioning import Toolset

ppath = 'analysis_driver.dataset.'
fake_proc = {'proc_id': 'a_proc_id', 'dataset_type': 'test', 'dataset_name': 'test'}

patched_patch = patch(ppath + 'rest_communication.patch_entry')
patched_post = patch(ppath + 'rest_communication.post_entry')
patched_pid = patch(ppath + 'os.getpid', return_value=1)
patched_update = patch(ppath + 'MostRecentProc.update_entity')
patched_initialise = patch(ppath + 'MostRecentProc.initialise_entity')
patched_finish = patch(ppath + 'MostRecentProc.finish')
patched_is_ready = patch(ppath + 'Dataset._is_ready', return_value=None)
patched_stages = patch(
    ppath + 'Dataset.running_stages',
    new_callable=PropertyMock(return_value=['this', 'that', 'other'])
)
patched_get_run = patch(
    ppath + 'clarity.get_run',
    return_value=Mock(udf={'Run Status': 'RunStarted'})
)


def patched_get_doc(content=None):
    return patch(ppath + 'rest_communication.get_document', return_value=content or fake_proc)


def patched_get_docs(content=None):
    return patch(ppath + 'rest_communication.get_documents', return_value=content or [fake_proc])


def patched_datetime(time='now'):
    return patch(ppath + 'now', return_value=time)


def patched_required_yield(y=1000000000):
    return patch(ppath + 'rest_communication.get_document', return_value={'required_yield': y})


class TestDataset(TestAnalysisDriver):
    base_dir = os.path.join(TestAnalysisDriver.assets_path, 'dataset_scanner')

    def setUp(self):
        self.setup_dataset()
        self.dataset._ntf = Mock()

    def test_dataset_status(self):
        assert self.dataset.most_recent_proc.entity.get('status') is None

        self.dataset.most_recent_proc.entity['status'] = 'a_status'
        assert self.dataset.dataset_status == 'a_status'

    @patched_get_docs([{'stage_name': 'this', 'date_started': 'now'}, {'stage_name': 'that', 'date_started': 'then'}])
    def test_stages(self, mocked_get):
        assert self.dataset.running_stages == ['this', 'that']
        mocked_get.assert_called_with(
            'analysis_driver_stages',
            all_pages=True,
            where={'analysis_driver_proc': 'a_proc_id', 'date_finished': None}
        )

    @patch(ppath + 'rest_communication.get_document', return_value={'run_id': 'a_run'})  # for RunDataset.initialise_entity
    @patched_initialise
    @patch(ppath + 'MostRecentProc.start')
    @patch(ppath + 'MostRecentProc.retrieve_entity')
    @patch(ppath + 'SampleDataset.lims_ntf', new_callable=PropertyMock)
    def test_start(self, mocked_lims_ntf, mocked_retrieve_entity, mocked_start, mocked_update, mocked_get_doc):
        self.dataset.most_recent_proc.entity['status'] = c.DATASET_PROCESSING
        with patched_get_run, patched_is_ready:
            with pytest.raises(AssertionError):
                self.dataset.start()
            del self.dataset.most_recent_proc.entity['status']
            self.dataset.start()
        for m in (mocked_update, self.dataset.ntf.start_pipeline, mocked_start):
            assert m.call_count == 1

    @patched_finish
    @patch(ppath + 'MostRecentProc.retrieve_entity')
    @patch(ppath + 'SampleDataset.lims_ntf', new_callable=PropertyMock)
    def test_succeed(self, mocked_lims_ntf, mocked_retrieve_entity, mocked_finish):
        with patched_get_run, patched_is_ready:
            with pytest.raises(AssertionError):
                self.dataset.succeed()
            self.dataset.most_recent_proc.entity['status'] = c.DATASET_PROCESSING
            self.dataset.succeed()
        self.dataset.ntf.end_pipeline.assert_called_with(0)
        mocked_finish.assert_called_with(c.DATASET_PROCESSED_SUCCESS)

    @patched_finish
    @patch(ppath + 'MostRecentProc.retrieve_entity')
    @patch(ppath + 'SampleDataset.lims_ntf', new_callable=PropertyMock)
    def test_fail(self, mocked_lims_ntf, mocked_retrieve_entity, mocked_finish):
        with patched_get_run, patched_is_ready:
            with pytest.raises(AssertionError):
                self.dataset.succeed()

        self.dataset.most_recent_proc.entity['status'] = c.DATASET_PROCESSING
        self.dataset.fail(1)
        self.dataset.ntf.end_pipeline.assert_called_with(1)
        mocked_finish.assert_called_with(c.DATASET_PROCESSED_FAIL)

    @patched_finish
    @patch(ppath + 'MostRecentProc.retrieve_entity')
    @patch(ppath + 'SampleDataset.lims_ntf', new_callable=PropertyMock)
    def test_abort(self, mocked_lims_ntf, mocked_retrieve_entity, mocked_finish):
        self.dataset.abort()
        mocked_finish.assert_called_with(c.DATASET_ABORTED)

    @patch(ppath + 'MostRecentProc.change_status')
    @patch(ppath + 'MostRecentProc.retrieve_entity')
    def test_reset(self, mocked_retrieve_entity, mocked_change_status):
        with patch(ppath + 'Dataset.terminate'):
            self.dataset.reset()
        mocked_change_status.assert_called_with(c.DATASET_REPROCESS)

    @patch(ppath + 'sleep')
    @patch(ppath + 'os.kill')
    @patch(ppath + 'Dataset._pid_running', side_effect=[True, False])
    def test_terminate(self, mocked_running, mocked_kill, mocked_sleep):
        self.dataset.most_recent_proc._entity['pid'] = 1337
        with patch(ppath + 'Dataset._pid_valid', return_value=False) as mocked_valid:
            self.dataset._terminate(signal.SIGUSR2)
            mocked_valid.assert_called_with(1337)
            assert all(m.call_count == 0 for m in (mocked_kill, mocked_running, mocked_sleep))

        with patch(ppath + 'Dataset._pid_valid', return_value=True) as mocked_valid:
            self.dataset.terminate()
            mocked_valid.assert_called_with(1337)
            mocked_kill.assert_called_with(1337, signal.SIGUSR2)
            assert mocked_running.call_count == 2
            assert mocked_sleep.call_count == 1

    def test_pid_valid(self):
        cmdlinefile = os.path.join(TestAnalysisDriver.assets_path, 'example.pid')
        if os.path.isfile(cmdlinefile):
            os.remove(cmdlinefile)

        with patch(ppath + 'os.path.join', return_value=cmdlinefile):
            assert self.dataset._pid_valid(1337) is False

            with open(cmdlinefile, 'w') as f:
                f.write(modules['__main__'].__file__ + '\n')

            assert self.dataset._pid_valid(1337) is True

    @patch(ppath + 'MostRecentProc.start_stage')
    def test_start_stage(self, mocked_start_stage):
        self.dataset.start_stage('a_stage')
        mocked_start_stage.assert_called_with('a_stage')
        self.dataset.ntf.start_stage.assert_called_with('a_stage')

    @patch(ppath + 'MostRecentProc.end_stage')
    def test_end_stage(self, mocked_end_stage):
        self.dataset.end_stage('a_stage')
        mocked_end_stage.assert_called_with('a_stage', 0)
        self.dataset.ntf.end_stage.assert_called_with('a_stage', 0)

    def test_report(self):
        with patched_stages:
            assert self.dataset.report() == 'test_dataset -- this, that, other'

    def setup_dataset(self):
        with patched_get_docs():
            self.dataset = Dataset(
                'test_dataset',
                {'proc_id': 'a_proc_id', 'date_started': 'now', 'dataset_name': 'None', 'dataset_type': 'None'}
            )

    def test_raise_exceptions(self):
        self.dataset.register_exception(Mock(stage_name='task1'), SequencingRunError('RunErrored'))
        # Give time for the exception to be pickled.
        # From the documentation:
        # 'After putting an object on an empty queue there may be an infinitesimal delay before the queue’s empty()
        # method returns False and get_nowait() can return without raising queue.Empty.'
        time.sleep(.1)
        with pytest.raises(SequencingRunError):
            self.dataset.raise_exceptions()

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
            'toolset_version': 1  # as per tests/assets/tool_versioning.yaml
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
    def test_is_ready(self):
        with patched_get_docs():
            d = RunDataset('dataset_ready')
            assert d._is_ready()

    def test_dataset_status(self):
        with patched_get_run:
            super().test_dataset_status()
            del self.dataset.most_recent_proc.entity['status']
            assert self.dataset.dataset_status == c.DATASET_READY

    def setup_dataset(self):
        self.dataset = RunDataset(
            'test_dataset',
            most_recent_proc={'proc_id': 'a_proc_id', 'date_started': 'now',
                              'dataset_name': 'None', 'dataset_type': 'None'}
        )

    def test_run_elements_from_lims(self):
        d = RunDataset('test_dataset')
        with patch('analysis_driver.dataset.clarity.get_run', return_value=MockedRunProcess(container=mocked_flowcell_non_pooling)), \
             patch.object(RunDataset, 'has_barcodes', new_callable=PropertyMock(return_value=False)):
            run_elements = d._run_elements_from_lims()
            assert len(set(r[c.ELEMENT_PROJECT_ID] for r in run_elements)) == 1
            assert len(set(r[c.ELEMENT_SAMPLE_INTERNAL_ID] for r in run_elements)) == 2
            barcodes_len = set(len(r[c.ELEMENT_BARCODE]) for r in run_elements)
            assert len(barcodes_len) == 1
            assert barcodes_len.pop() == 0

        d = RunDataset('test_dataset')
        with patch('egcg_core.clarity.get_run', return_value=MockedRunProcess(container=mocked_flowcell_pooling)), \
             patch.object(RunDataset, 'has_barcodes', new_callable=PropertyMock(return_value=True)):
            run_elements = d._run_elements_from_lims()
            assert len(set(r[c.ELEMENT_PROJECT_ID] for r in run_elements)) == 1
            assert len(set(r[c.ELEMENT_SAMPLE_INTERNAL_ID] for r in run_elements)) == 4
            barcodes_len = set(len(r[c.ELEMENT_BARCODE]) for r in run_elements)
            assert len(barcodes_len) == 1
            assert barcodes_len.pop() == 8

        d = RunDataset('test_dataset')
        with patch('analysis_driver.dataset.clarity.get_run', return_value=MockedRunProcess(container=mocked_flowcell_user_prepared)), \
             patch.object(RunDataset, 'has_barcodes', new_callable=PropertyMock(return_value=False)):
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



class TestSampleDataset(TestDataset):
    def test_dataset_status(self):
        super().test_dataset_status()
        del self.dataset.most_recent_proc.entity['status']
        assert self.dataset.dataset_status == c.DATASET_READY

    @patch(ppath + 'MostRecentProc.change_status')
    def test_force(self, mocked_change_status):
        self.dataset.force()
        mocked_change_status.assert_called_with(c.DATASET_FORCE_READY)

    def test_amount_data(self):
        assert self.dataset._amount_data() == 1500000000

    def test_data_threshold(self):
        self.dataset._data_threshold = None
        assert self.dataset.data_threshold == 1000000000

    def test_no_data_threshold(self):
        self.dataset._data_threshold = None
        with pytest.raises(AnalysisDriverError) as e:
            self.dataset._sample.pop('required_yield')
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
            'toolset_version': 1
        }

        with patched_get_doc(None), patched_patch as mocked_patch:
            assert self.dataset._pipeline_instruction() == exp
            mocked_patch.assert_called_with('projects', {'sample_pipeline': exp}, 'project_id', 'a_project')

    def test_is_ready(self):
        self.dataset._data_threshold = 20000000000
        assert not self.dataset._is_ready()
        self.dataset._data_threshold = 100
        assert self.dataset._is_ready()

    def test_report(self):
        expected_str = 'test_dataset -- this, that, other  (1500000000 / 1000000000  from a_run_id, another_run_id) (non useable run elements in a_run_id, another_run_id)'
        self.dataset._data_threshold = None
        with patched_required_yield(), patched_stages:
            assert self.dataset.report() == expected_str

    def setup_dataset(self):
        with patched_get_docs():
            self.dataset = SampleDataset(
                'test_dataset',
                most_recent_proc={'proc_id': 'a_proc_id', 'date_started': 'now',
                                  'dataset_name': 'None', 'dataset_type': 'None'}
            )
        self.dataset._run_elements = self.dataset._non_useable_run_elements = [
            {'run_element_id': 'run_element1', 'run_id': 'a_run_id', 'clean_q30_bases_r1': 120, 'clean_q30_bases_r2': 115, 'q30_bases_r1': 150, 'q30_bases_r2': 130, 'bases_r1': 200, 'bases_r2': 190},
            {'run_element_id': 'run_element2', 'run_id': 'another_run_id', 'clean_q30_bases_r1': 110, 'clean_q30_bases_r2': 135, 'q30_bases_r1': 170, 'q30_bases_r2': 150, 'bases_r1': 210, 'bases_r2': 205}
        ]
        self.dataset._sample = {
            'aggregated': {'clean_yield_in_gb': 1.5, 'run_ids': ['a_run_id', 'another_run_id'], 'clean_pc_q30': 85},
            'required_yield': 1000000000
        }
        self.dataset._data_threshold = 1000000000

    @patched_initialise
    @patch(ppath + 'MostRecentProc.start')
    @patch(ppath + 'MostRecentProc.retrieve_entity')
    @patch(ppath + 'SampleDataset.lims_ntf', new_callable=PropertyMock)
    def test_start(self, mocked_lims_ntf, mocked_retrieve_entity, mocked_start, mocked_update):
        with patched_get_run:
            self.dataset.start()
            mocked_lims_ntf.assert_has_calls([call(), call().start_sample_pipeline()])

    @patched_finish
    @patch(ppath + 'MostRecentProc.retrieve_entity')
    @patch(ppath + 'SampleDataset.lims_ntf', new_callable=PropertyMock)
    def test_succeed(self, mocked_lims_ntf, mocked_retrieve_entity, mocked_finish):
        with patched_get_run:
            self.dataset.most_recent_proc.entity['status'] = c.DATASET_PROCESSING
            self.dataset.succeed()
            mocked_lims_ntf.assert_has_calls([call(), call().assign_next_and_advance_step()])

    @patched_finish
    @patch(ppath + 'MostRecentProc.retrieve_entity')
    @patch(ppath + 'SampleDataset.lims_ntf', new_callable=PropertyMock)
    def test_fail(self, mocked_lims_ntf, mocked_retrieve_entity, mocked_finish):
        self.dataset.most_recent_proc.entity['status'] = c.DATASET_PROCESSING
        self.dataset.fail(1)
        self.assertTrue(mocked_lims_ntf.called)
        mocked_lims_ntf.assert_has_calls([call(), call().remove_sample_from_workflow()])

    def test_data_source(self):
        data_source = self.dataset.data_source
        assert data_source == ['run_element1', 'run_element2']


class TestMostRecentProc(TestAnalysisDriver):
    def setUp(self):
        with patched_get_docs():
            self.proc = MostRecentProc(
                NamedMock(type='test', endpoint='tests', real_name='test', pipeline_type='some kind of variant calling', data_source=None, id_field='test_id')
            )
        self.proc._entity = fake_proc.copy()
        self.proc.proc_id = 'a_proc_id'

    def test_rest_entity_pre_existing(self):
        assert self.proc.entity == fake_proc

    @patched_get_docs()
    @patched_initialise
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
        sample_proc = MostRecentProc(NamedMock(type='sample', real_name='sample1', pipeline_type='some kind of variant calling', data_source=['run_element1', 'run_element2']))
        sample_proc._entity = None
        entity = sample_proc.entity
        with patched_datetime():

            sample_proc.initialise_entity()
            mocked_post.assert_called_with('analysis_driver_procs',
                                           {'dataset_type': 'sample',
                                            'dataset_name': 'sample1',
                                            'data_source': ['run_element1', 'run_element2'],
                                            'proc_id': 'sample_sample1_now'})

        project_proc = MostRecentProc(NamedMock(type='project', real_name='project1', pipeline_type='project pipeline', data_source=['sample_1', 'sample_2']))
        project_proc._entity = None
        entity = project_proc.entity
        with patched_datetime():
            project_proc.initialise_entity()
            mocked_post.assert_called_with('analysis_driver_procs',
                                           {'data_source': ['sample_1', 'sample_2'],
                                            'dataset_name': 'project1',
                                            'proc_id': 'project_project1_now',
                                            'dataset_type': 'project'})

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
        with patched_pid:
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
