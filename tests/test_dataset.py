import os
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


class TestDataset(TestAnalysisDriver):
    base_dir = os.path.join(TestAnalysisDriver.assets_path, 'dataset_scanner')

    def setUp(self):
        self.dataset = Dataset('test_dataset')
        self.dataset._ntf = Mock()
        self.dataset.most_recent_proc = Mock()

    def test_dataset_status(self):
        self.dataset.most_recent_proc = {'status': 'a_status'}
        assert self.dataset.dataset_status == 'a_status'

        self.dataset.most_recent_proc = {}
        with patch.object(self.dataset.__class__, '_is_ready', return_value=True) as mocked_ready:
            assert self.dataset.dataset_status == c.DATASET_READY
            mocked_ready.return_value = False
            assert self.dataset.dataset_status == c.DATASET_NEW

    @patched_get_docs([{'stage_name': 'this', 'date_started': 'now'}, {'stage_name': 'that', 'date_started': 'then'}])
    def test_running_stages(self, mocked_get):
        self.dataset.most_recent_proc = {'proc_id': 'a_proc_id'}
        assert self.dataset.running_stages == ['this', 'that']
        mocked_get.assert_called_with(
            'analysis_driver_stages',
            all_pages=True,
            where={'analysis_driver_proc': 'a_proc_id', 'date_finished': None}
        )

    def test_start(self):
        self.dataset.most_recent_proc.get.return_value = c.DATASET_RESUME
        self.dataset.start()
        assert self.dataset.most_recent_proc.retrieve_entity.call_count == 2  # 1 for assert_status, 1 for start
        assert self.dataset.most_recent_proc.start.call_count == 1
        assert self.dataset.ntf.start_pipeline.call_count == 1

        self.dataset.most_recent_proc.get.return_value = c.DATASET_READY
        self.dataset.start()
        assert self.dataset.most_recent_proc.initialise_entity.call_count == 1

    def test_succeed(self):
        self.dataset.most_recent_proc.get.return_value = c.DATASET_PROCESSING
        self.dataset.succeed()
        self.dataset.most_recent_proc.finish.assert_called_with(c.DATASET_PROCESSED_SUCCESS)
        self.dataset.ntf.end_pipeline.assert_called_with(0)
        self.dataset.most_recent_proc.get.return_value = {'status': c.DATASET_READY}
        with self.assertRaises(AssertionError):
            self.dataset.succeed()

    def test_fail(self):
        self.dataset.most_recent_proc.get.return_value = c.DATASET_PROCESSING
        self.dataset.fail(1)
        self.dataset.most_recent_proc.finish.assert_called_with(c.DATASET_PROCESSED_FAIL)
        self.dataset.ntf.end_pipeline.assert_called_with(1)

    def test_abort(self):
        self.dataset.abort()
        self.dataset.most_recent_proc.finish.assert_called_with(c.DATASET_ABORTED)

    def test_resume(self):
        with patch.object(self.dataset.__class__, 'terminate'):
            self.dataset.resume()
        self.dataset.most_recent_proc.change_status.assert_called_with(c.DATASET_RESUME)

    def test_reset(self):
        with patch.object(self.dataset.__class__, 'terminate'):
            self.dataset.reset()
        self.dataset.most_recent_proc.change_status.assert_called_with(c.DATASET_REPROCESS)

    def test_force(self):
        self.dataset.force()
        self.dataset.most_recent_proc.change_status.assert_called_with(c.DATASET_FORCE_READY)

    @patch(ppath + 'sleep')
    @patch('os.kill')
    def test_terminate(self, mocked_kill, mocked_sleep):
        # no pid - nothing happens
        self.dataset.most_recent_proc = {'pid': None}
        self.dataset._terminate(1)
        assert mocked_kill.call_count == 0

        # pid exists in database but has been finished or recycled - nothing happens
        self.dataset.most_recent_proc = {'pid': 1337}
        with patch.object(self.dataset.__class__, '_pid_valid', return_value=False) as mocked_valid:
            self.dataset._terminate(1)
            mocked_valid.assert_called_with(1337)
            assert mocked_kill.call_count == 0

            mocked_valid.return_value = True
            with patch.object(self.dataset.__class__, '_pid_running', side_effect=[True, False]):
                self.dataset._terminate(1)

        mocked_kill.assert_called_with(1337, 1)
        mocked_sleep.assert_called_with(1)

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

    def test_report(self):
        self.dataset.most_recent_proc = {'pid': 1337}
        with patched_stages:
            assert self.dataset.report() == 'test_dataset (1337) -- this, that, other'

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
        self.dataset._processing_instruction = Mock(
            return_value={'toolset_type': 'a_toolset', 'toolset_version': 3, 'name': 'qc'}
        )
        with patch(ppath + 'toolset') as mocked_toolset:
            self.dataset.resolve_pipeline_and_toolset()

        mocked_toolset.configure.assert_called_with(
            'a_toolset',
            3,
            'tests/assets/jobs/test_dataset/program_versions.yaml'
        )

    @patch(ppath + 'toolset', new=Toolset())
    def test_processing_instruction(self):
        pipeline = NamedMock(real_name='qc', toolset_type='non_human_sample_processing')
        with patch.object(self.dataset.__class__, 'pipeline', new=pipeline):
            assert self.dataset._processing_instruction() == {
                'name': 'qc',
                'toolset_type': 'non_human_sample_processing',
                'toolset_version': 1  # as per tests/assets/tool_versioning.yaml
            }


mocked_lane_user_prep_artifact1 = NamedMock(real_name='art1', reagent_labels=[], samples=[MockedSample(real_name='sample1')])
mocked_lane_user_prep_artifact2 = NamedMock(real_name='art2', reagent_labels=[], samples=[MockedSample(real_name='sample2')])

mocked_lane_artifact1 = NamedMock(real_name='art1', reagent_labels=['D701-D502 (ATTACTCG-ATAGAGGC)'], samples=[MockedSample(real_name='sample1', udf={})])
mocked_lane_artifact2 = NamedMock(real_name='art2', reagent_labels=['D702-D502 (TCCGGAGA-ATAGAGGC)'], samples=[MockedSample(real_name='sample2', udf={'Rapid Analysis': 'Yes'})])
mocked_lane_artifact3 = NamedMock(real_name='art3', reagent_labels=['D703-D502 (CGCTCATT-ATAGAGGC)'], samples=[MockedSample(real_name='sample3', udf={})])
mocked_lane_artifact4 = NamedMock(real_name='art4', reagent_labels=['D704-D502 (GAGATTCC-ATAGAGGC)'], samples=[MockedSample(real_name='sample4', udf={})])
mocked_lane_artifact5 = NamedMock(real_name='art5', reagent_labels=['D705-D502 (ATTCAGAA-ATAGAGGC)'], samples=[MockedSample(real_name='sample5', udf={})])
mocked_lane_artifact6 = NamedMock(real_name='art6', reagent_labels=['D706-D502 (GAATTCGT-ATAGAGGC)'], samples=[MockedSample(real_name='sample6', udf={})])
mocked_lane_artifact7 = NamedMock(real_name='art7', reagent_labels=['D706-D502 (GAATTCGA-ATAGAGGC)'], samples=[MockedSample(real_name='sample7', udf={})])
mocked_lane_artifact8 = NamedMock(real_name='art8', reagent_labels=['D706-D502 (GAATTCGG-ATAGAGGC)'], samples=[MockedSample(real_name='sample8', udf={})])
mocked_idt_artifact = NamedMock(
    real_name='idt_art',
    reagent_labels=['001A IDT-ILMN TruSeq DNA-RNA UD 96 Indexes  Plate_UDI0001 (CCGCGGTT-AGCGCTAG)'],
    samples=[MockedSample(real_name='idt_sample')]
)
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
    '3:1': mocked_lane_artifact3,
    '4:1': mocked_lane_artifact4,
    '5:1': mocked_lane_artifact5,
    '6:1': mocked_lane_artifact6,
    '7:1': mocked_lane_artifact7,
    '8:1': mocked_lane_artifact8
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

mocked_flowcell_idt = Mock(placements={'%s:1' % i: mocked_idt_artifact for i in range(1, 9)})


class TestRunDataset(TestDataset):
    def setUp(self):
        self.dataset = RunDataset('test_dataset')
        self.dataset._ntf = Mock()
        self.dataset.most_recent_proc = Mock()

    def test_start(self):
        with patched_get_doc({'some': 'content'}):
            super().test_start()

    def test_is_ready(self):
        with patched_get_docs():
            d = RunDataset('dataset_ready')
            assert d._is_ready()

    @patch('analysis_driver.dataset.clarity.get_run')
    def test_rapid_samples_by_lane(self, mocked_get_run):
        mocked_get_run.return_value = MockedRunProcess(container=mocked_flowcell_pooling)
        assert self.dataset.rapid_samples_by_lane == {}

        self.dataset._rapid_samples_by_lane = None
        mocked_get_run.return_value.container = mocked_flowcell_non_pooling
        assert self.dataset.rapid_samples_by_lane == {'2': {'sample_id': 'sample2', 'Rapid Analysis': 'Yes', 'project_id': '10015AT'}}

    def test_run_elements_from_lims(self):
        def patched_run(fake_flowcell):
            return patch(
                'analysis_driver.dataset.clarity.get_run',
                return_value=MockedRunProcess(container=fake_flowcell)
            )

        def patched_barcodes(has_barcodes):
            return patch.object(RunDataset, 'has_barcodes', new_callable=PropertyMock(return_value=has_barcodes))

        d = RunDataset('test_dataset')
        with patched_run(mocked_flowcell_non_pooling), patched_barcodes(False):
            run_elements = d._run_elements_from_lims()
            assert len(set(r[c.ELEMENT_PROJECT_ID] for r in run_elements)) == 1
            assert len(set(r[c.ELEMENT_SAMPLE_INTERNAL_ID] for r in run_elements)) == 8
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

        d = RunDataset('test_dataset')
        with patched_run(mocked_flowcell_idt), patched_barcodes(True):
            run_elements = d._run_elements_from_lims()
            assert len(set(r[c.ELEMENT_PROJECT_ID] for r in run_elements)) == 1
            assert len(set(r[c.ELEMENT_SAMPLE_INTERNAL_ID] for r in run_elements)) == 1
            barcodes_len = set(len(r[c.ELEMENT_BARCODE]) for r in run_elements)
            assert len(barcodes_len) == 1
            assert barcodes_len.pop() == 8

    @patch('builtins.open')
    def test_generate_samplesheet(self, mocked_open):
        fake_run_elements = [
            {'lane': '2', 'sample_id': 'sample_2', 'library_id': 'lib_2', 'project_id': 'p', 'barcode': 'CTGA'},
            {'lane': '1', 'sample_id': 'sample_1', 'library_id': 'lib_1', 'project_id': 'p', 'barcode': 'ATGC'}
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

    def test_sample_sheet_file(self):
        self.dataset.input_dir = os.path.join(self.assets_path)
        sample_sheet_file = os.path.join(self.dataset.input_dir, 'SampleSheet_analysis_driver.csv')
        with patch.object(RunDataset, '_generate_samplesheet') as mgenerate:
            _ = self.dataset.sample_sheet_file
        mgenerate.assert_called_once_with(sample_sheet_file)

    def test_sample_sheet_file_exists(self):
        self.dataset.input_dir = os.path.join(self.assets_path)
        sample_sheet_file = os.path.join(self.dataset.input_dir, 'SampleSheet_analysis_driver.csv')
        open(sample_sheet_file, 'w').close()
        with patch.object(RunDataset, '_generate_samplesheet') as mgenerate:
            _ = self.dataset.sample_sheet_file
        assert mgenerate.call_count == 0
        os.remove(sample_sheet_file)


class TestSampleDataset(TestDataset):
    def setUp(self):
        self.dataset = SampleDataset('test_dataset')
        self.dataset._ntf = Mock()
        self.dataset._lims_ntf = Mock()
        self.dataset.most_recent_proc = Mock()

        self.dataset._run_elements = self.dataset._non_useable_run_elements = [
            {'run_element_id': 'run_element1', 'run_id': 'a_run_id', 'clean_q30_bases_r1': 120,
             'clean_q30_bases_r2': 115, 'q30_bases_r1': 150, 'q30_bases_r2': 130, 'bases_r1': 200, 'bases_r2': 190},
            {'run_element_id': 'run_element2', 'run_id': 'another_run_id', 'clean_q30_bases_r1': 110,
             'clean_q30_bases_r2': 135, 'q30_bases_r1': 170, 'q30_bases_r2': 150, 'bases_r1': 210, 'bases_r2': 205}
        ]
        self.dataset._sample = {
            'aggregated': {'clean_yield_in_gb': 1.5, 'run_ids': ['a_run_id', 'another_run_id'], 'clean_pc_q30': 85,
                           'from_run_elements': {'mean_coverage': 32}},
            'required_yield': 1000000000,
            'required_coverage': 30,
            'project_id': 'a_project'
        }
        self.dataset._species = 'Teleogryllus oceanicus'
        self.dataset._genome_version = '1'
        self.dataset._data_threshold = 1000000000

    @patched_get_doc({'project_whitelist': ['another_project']})
    def test_genome_dict(self, mocked_get):
        with self.assertRaises(AnalysisDriverError) as e:
            _ = self.dataset.genome_dict

        assert str(e.exception) == 'Project ID a_project not in whitelist for reference genome 1'

        mocked_get.return_value = {
            'project_whitelist': ['a_project'],
            'data_files': {'fasta': 'Teleogryllus_oceanicus/v1/ref.fa'}
        }
        assert self.dataset.genome_dict == {
            'project_whitelist': ['a_project'],
            'data_files': {'fasta': 'path/to/genomes_dir/Teleogryllus_oceanicus/v1/ref.fa'}
        }

    def test_amount_data(self):
        assert self.dataset._amount_data() == 1500000000

    def test_required_yield_threshold(self):
        assert self.dataset.required_yield_threshold == 1000000000

    def test_no_required_yield_threshold(self):
        with pytest.raises(AnalysisDriverError) as e:
            self.dataset._sample.pop('required_yield')
            _ = self.dataset.required_yield_threshold
        assert 'Could not find required yield threshold' in str(e)

    def test_required_coverage_threshold(self):
        assert self.dataset.required_coverage_threshold == 30

    def test_no_required_coverage_threshold(self):
        with pytest.raises(AnalysisDriverError) as e:
            self.dataset._sample.pop('required_coverage')
            _ = self.dataset.required_coverage_threshold
        assert 'Could not find required coverage threshold' in str(e)

    @patched_patch
    @patch(ppath + 'toolset', new=Toolset())
    @patched_get_doc({'sample_pipeline': 'some_data'})
    @patch.object(SampleDataset, 'pipeline', new=NamedMock(real_name='qc', toolset_type='non_human_sample_processing'))
    def test_processing_instruction(self, mocked_get, mocked_patch):
        assert self.dataset._processing_instruction() == 'some_data'
        assert mocked_patch.call_count == 0

        exp = {
            'name': 'qc',
            'toolset_type': 'non_human_sample_processing',
            'toolset_version': 1  # as per tests/assets/tool_versioning.yaml
        }
        mocked_get.return_value = {}
        assert self.dataset._processing_instruction() == exp
        mocked_patch.assert_called_with(
            'projects',
            {'sample_pipeline': {'name': 'qc', 'toolset_type': 'non_human_sample_processing', 'toolset_version': 1}},
            'project_id',
            'a_project'
        )

    def test_is_ready(self):
        self.dataset.sample['required_yield'] = 20000000000
        self.dataset.sample['required_coverage'] = 60
        assert not self.dataset._is_ready()
        self.dataset.sample['required_yield'] = 100
        assert self.dataset._is_ready()

    def test_report(self):
        self.dataset.most_recent_proc = {'pid': 1337}
        with patched_stages:
            assert self.dataset.report() == (
                'test_dataset (1337) -- this, that, other  (yield: 1500000000 / 1000000000 and coverage: 32 / 30 '
                'from a_run_id, another_run_id) (non useable run elements in a_run_id, another_run_id)'
            )

    def test_start(self):
        super().test_start()
        assert self.dataset.lims_ntf.start_sample_pipeline.called

    def test_succeed(self):
        super().test_succeed()
        assert self.dataset.lims_ntf.assign_next_and_advance_step.called

    def test_fail(self):
        self.dataset.most_recent_proc.get.return_value = c.DATASET_PROCESSING
        self.dataset.fail(1)
        assert self.dataset.lims_ntf.remove_sample_from_workflow.called

    def test_data_source(self):
        data_source = self.dataset.data_source
        assert data_source == ['run_element1', 'run_element2']


class TestProjectDataset(TestDataset):
    def setUp(self):
        self.dataset = ProjectDataset('test_dataset')
        self.dataset._ntf = Mock()
        self.dataset.most_recent_proc = Mock()

    def test_samples_processed(self):
        with patched_get_docs(['sample1', 'sample2']) as mgetdoc:
            assert self.dataset.samples_processed == ['sample1', 'sample2']
            mgetdoc.assert_called_with(
                'samples', all_pages=True,
                where={'project_id': 'test_dataset', 'aggregated.most_recent_proc.status': 'finished'}
            )

    def test_get_processed_gvcfs(self):
        self.dataset._samples_processed = [
            {
                'sample_id': 'sample_1',
                'user_sample_id': 'uid_1',
                'aggregated': {'most_recent_proc': {'pipeline_used': {'name': 'bcbio'}}}
            },
            {
                'sample_id': 'sample_2',
                'user_sample_id': 'uid_2',
                'aggregated': {'most_recent_proc': {'pipeline_used': {'name': 'not_bcbio'}}}
            }
        ]

        with patch(ppath + 'find_file', new=lambda *args: '/'.join(args)):
            obs = self.dataset.get_processed_gvcfs()
            assert obs == [
                'tests/assets/test_projects/test_dataset/sample_1/uid_1.g.vcf.gz'
            ]

    @patch(ppath + 'clarity.get_project', return_value=None)
    def test_number_of_samples(self, mocked_project):
        assert self.dataset.number_of_samples == -1  # no project data

        self.dataset._number_of_samples = None
        mocked_project.return_value = Mock(udf={})
        assert self.dataset.number_of_samples == -1  # quoted samples missing

        self.dataset._number_of_samples = None
        mocked_project.return_value = Mock(udf={'Number of Quoted Samples': 3})
        assert self.dataset.number_of_samples == 3

    @patch(ppath + 'clarity.get_species_from_sample', return_value='species_from_clarity')
    def test_species(self, mocked_species):
        self.dataset._samples_processed = [
            {'sample_id': 'sample_1', 'species_name': 'this'},
            {'sample_id': 'sample_2'}
        ]

        # species = 'this', 'species_from_clarity'
        with self.assertRaises(AnalysisDriverError):
            _ = self.dataset.species

        # species = 'this', 'this'
        self.dataset._species = None
        self.dataset._samples_processed[1]['species_name'] = 'this'
        assert self.dataset.species == 'this'

    @patched_get_doc({'default_version': '1'})
    @patch(ppath + 'clarity.get_sample', return_value=Mock(udf={'Genome Version': '2'}))
    def test_genome_version(self, mocked_sample, mocked_get):
        self.dataset._samples_processed = [{'sample_id': 'sample_1'}]
        self.dataset._species = 'this'
        assert self.dataset.genome_version == '2'
        mocked_sample.assert_called_with('sample_1')
        assert mocked_get.call_count == 0

        self.dataset._genome_version = None
        self.dataset._reference_genome = None
        mocked_sample.return_value = Mock(udf={})
        assert self.dataset.genome_version == '1'
        mocked_get.assert_called_with('species', where={'name': 'this'})


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
        self.proc.dataset.type = 'sample'
        self.proc.dataset.genome_version = 'Tthi_1.0'
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
                },
                'genome_used': 'Tthi_1.0'
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
