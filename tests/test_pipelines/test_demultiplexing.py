import os
import shutil
import pytest
from os.path import join
from unittest.mock import Mock, patch
from egcg_core.constants import ELEMENT_PROJECT_ID, ELEMENT_LANE, ELEMENT_SAMPLE_INTERNAL_ID
from tests.test_analysisdriver import TestAnalysisDriver, NamedMock
from analysis_driver.exceptions import SequencingRunError
from analysis_driver.util import bash_commands
from analysis_driver.pipelines import demultiplexing as dm

patch_executor = patch(
    'analysis_driver.pipelines.demultiplexing.executor.execute',
    return_value=Mock(join=Mock(return_value=0))
)

patched_ref_genome = patch('analysis_driver.dataset.RunDataset.reference_genome', return_value='a_genome.fa')


class TestPhixDetection(TestAnalysisDriver):
    def test_run(self):
        dataset = NamedMock(real_name='test')
        f = dm.PhixDetection(dataset=dataset)
        fake_fastq_pairs = [
            [('L1_R1_001.fastq.gz', 'L1_R2_001.fastq.gz')], [('L2_R1_001.fastq.gz', 'L2_R2_001.fastq.gz')]
        ]
        fake_genome_response = {
            "_updated": "30_11_2018_15:13:43",
            "assembly_name": "phix174",
            "analyses_supported": ["qc"],
            "data_source": "",
            "_links": {"self": {"title": "genome", "href": "genomes/phix174"}},
            "_etag": "175b41e3909a93a8298ac1d5d4dfc7292df4b580",
            "data_files": {"fasta": "path/to/phix.fa"},
            "_created": "30_11_2018_15:13:43",
            "species": "PhiX",
            "genome_size": 5386,
            "_id": "5c0153a716a5772f9e9cfdcc"
        }
        patch_run_crawler = patch('analysis_driver.pipelines.demultiplexing.RunCrawler', autospec=True)
        patch_find = patch('egcg_core.util.find_all_fastq_pairs', side_effect=fake_fastq_pairs)
        patch_get_document = patch('egcg_core.rest_communication.get_document',
                                   return_value=fake_genome_response)
        with patch_find, patch_executor as pexecute, patch_run_crawler as prun_crawler, patch_get_document:
            assert f._run() == 0

            assert pexecute.call_args[0][0] == (
                'set -o pipefail; path/to/bwa_1.1 mem -t 16 path/to/genomes_dir/path/to/phix.fa L1_R1_001.fastq.gz | '
                'path/to/samtools_1.3.1 view -F 4 | cut -f 1 | sort -u > L1_phix_read_name.list'
            )
            prun_crawler.assert_called_with(
                dataset, run_dir='tests/assets/jobs/test/fastq', stage=prun_crawler.STAGE_CONVERSION
            )


class TestFastqFilter(TestAnalysisDriver):
    def test_run(self):
        run_info = Mock(
            reads=Mock(
                upstream_read=Mock(attrib={'NumCycles': '151'}),
                downstream_read=Mock(attrib={'NumCycles': '151'}),
                index_lengths=[8]
            )
        )
        dataset = NamedMock(
            real_name='test',
            lane_metrics=[
                {'lane_number': 3, 'aggregated': {'pc_q30': 73}},
                {'lane_number': 4, 'aggregated': {'pc_q30': 73}}
            ],
            run_info=run_info
        )
        f = dm.FastqFilter(dataset=dataset)

        fake_fastq_pairs = [
            [('L1_R1_001.fastq.gz', 'L1_R2_001.fastq.gz')], [('L2_R1_001.fastq.gz', 'L2_R2_001.fastq.gz')],
            [('L3_R1_001.fastq.gz', 'L3_R2_001.fastq.gz')], [('L4_R1_001.fastq.gz', 'L4_R2_001.fastq.gz')],
            [('L5_R1_001.fastq.gz', 'L5_R2_001.fastq.gz')], [('L6_R1_001.fastq.gz', 'L6_R2_001.fastq.gz')],
            [('L7_R1_001.fastq.gz', 'L7_R2_001.fastq.gz')], [('L8_R1_001.fastq.gz', 'L8_R2_001.fastq.gz')]
        ]
        patch_detector = patch('analysis_driver.pipelines.demultiplexing.BadTileCycleDetector')
        patch_find = patch('analysis_driver.pipelines.demultiplexing.find_all_fastq_pairs_for_lane',
                           side_effect=fake_fastq_pairs + fake_fastq_pairs)
        with patch_find, patch_executor as pexecute, patch_detector as pdetector:
            instance = pdetector.return_value
            instance.detect_bad_tiles.return_value = {3: [1101]}
            instance.detect_bad_cycles.return_value = {4: [310, 308, 307, 309]}
            assert f._run() == 0

            expected_call_l2 = (
                'run_filterer in_place L2_R1_001.fastq.gz L2_R2_001.fastq.gz L2_R1_001_filtered.fastq.gz '
                'L2_R2_001_filtered.fastq.gz L2_R1_001_filtered.fastq L2_R2_001_filtered.fastq L2_fastqfilterer.stats '
                'L2_phix_read_name.list L2_R1_001.fastq_discarded L2_R2_001.fastq_discarded'
            )
            expected_call_l3 = (
                'run_filterer keep_originals L3_R1_001.fastq.gz L3_R2_001.fastq.gz L3_R1_001_filtered.fastq.gz '
                'L3_R2_001_filtered.fastq.gz L3_R1_001_filtered.fastq L3_R2_001_filtered.fastq L3_fastqfilterer.stats '
                'L3_phix_read_name.list L3_R1_001.fastq_discarded L3_R2_001.fastq_discarded --remove_tiles 1101'
            )
            expected_call_l4 = (
                'run_filterer keep_originals L4_R1_001.fastq.gz L4_R2_001.fastq.gz L4_R1_001_filtered.fastq.gz '
                'L4_R2_001_filtered.fastq.gz L4_R1_001_filtered.fastq L4_R2_001_filtered.fastq L4_fastqfilterer.stats '
                'L4_phix_read_name.list L4_R1_001.fastq_discarded L4_R2_001.fastq_discarded --trim_r2 147'
            )
            assert expected_call_l2 == pexecute.call_args[0][1]
            assert expected_call_l3 == pexecute.call_args[0][2]
            assert expected_call_l4 == pexecute.call_args[0][3]


class TestWaitForRead2(TestAnalysisDriver):
    ppath = 'analysis_driver.quality_control.interop_metrics.get_last_cycles_with_existing_bcls'

    @patch('time.sleep')
    def test_run(self, mocked_sleep):
        # Run info states 150 cycles for first read, 8 index cycles, and 50 cycles for second read = 208
        run_info = Mock(reads=Mock(upstream_read=Mock(attrib={'NumCycles': '150'}), index_lengths=[8]))
        dataset = NamedMock(real_name='testrun', run_info=run_info, input_dir='path/to/input',
                            lims_run=Mock(udf={'Run Status': 'RunStarted'}))

        self.stage = dm.WaitForRead2(dataset=dataset)

        with patch(self.ppath, return_value=310) as mcycles:
            assert mocked_sleep.call_count == 0
            assert self.stage._run() == 0
            assert mcycles.call_count == 1
            mcycles.assert_called_once_with('path/to/input')

        # get_last_cycles_with_existing_bcls states first 207, then 208 cycles done
        with patch(self.ppath, side_effect=[207, 208]) as mcycles:
            assert self.stage._run() == 0
            assert mcycles.call_count == 2
            mcycles.assert_called_with('path/to/input')
            mocked_sleep.assert_called_with(1200)

    def test_run_aborted(self):
        # Run info states 150 cycles for first read, 8 index cycles, and 50 cycles for second read = 208
        run_info = Mock(reads=Mock(upstream_read=Mock(attrib={'NumCycles': '150'}), index_lengths=[8]))
        dataset = NamedMock(real_name='testrun', run_info=run_info, input_dir='path/to/input',
                            lims_run=Mock(udf={'Run Status': 'RunAborted'}))

        # get_last_cycles_with_existing_bcls states 208 cycles done
        with patch(self.ppath, return_value=208) as mcycles:
            with pytest.raises(SequencingRunError):
                self.stage = dm.WaitForRead2(dataset=dataset)
                self.stage._run()
            assert mcycles.call_count == 1
            mcycles.assert_called_once_with('path/to/input')


class TestPostDemultiplexing(TestAnalysisDriver):
    @staticmethod
    def _touch(input_file):
        open(input_file, 'w').close()

    def setup_stage(self, stage_class):
        self.run_elements = [
            {ELEMENT_PROJECT_ID: 'project1', ELEMENT_SAMPLE_INTERNAL_ID: 'testsample1', ELEMENT_LANE: 1},
            {ELEMENT_PROJECT_ID: 'project1', ELEMENT_SAMPLE_INTERNAL_ID: 'testsample2', ELEMENT_LANE: 2}
        ]
        dataset = NamedMock(
            real_name='testrun',
            run_elements=self.run_elements
        )
        self.stage = stage_class(dataset=dataset)
        os.makedirs(self.stage.fastq_dir, exist_ok=True)
        os.makedirs(self.stage.fastq_intermediate_dir, exist_ok=True)

        for run_element in self.run_elements:
            s_dir = join(self.stage.fastq_intermediate_dir, run_element[ELEMENT_PROJECT_ID], run_element[ELEMENT_SAMPLE_INTERNAL_ID])
            os.makedirs(s_dir, exist_ok=True)
            self._touch(join(s_dir, 'someid_L00%s_R1_001.fastq.gz' % run_element[ELEMENT_LANE]))
            self._touch(join(s_dir, 'someid_L00%s_R2_001.fastq.gz' % run_element[ELEMENT_LANE]))

    def setup_post_alignment(self):
        for run_element in self.run_elements:
            s_dir = join(self.stage.job_dir, 'alignment', run_element[ELEMENT_PROJECT_ID], run_element[ELEMENT_SAMPLE_INTERNAL_ID])
            os.makedirs(s_dir, exist_ok=True)
            self._touch(join(s_dir, 'someid_L00%s.bam' % run_element[ELEMENT_LANE]))

    def tearDown(self):
        shutil.rmtree(self.stage.job_dir)


class TestBwaAlignMulti(TestPostDemultiplexing):
    def setUp(self):
        self.setup_stage(dm.BwaAlignMulti)

    def test_run(self):
        self.stage.dataset.reference_genome.return_value = 'a_genome.fa'
        with patch_executor as mock_executor:
            assert self.stage._run() == 0
            cmds = []
            for run_element in self.run_elements:
                cmds.append(bash_commands.bwa_mem_biobambam(
                    self.stage.fastq_pair(run_element),
                    'a_genome.fa',
                    self.stage.bam_path(run_element),
                    {'ID': '1', 'SM': run_element.get(ELEMENT_SAMPLE_INTERNAL_ID), 'PL': 'illumina'},
                    thread=6
                ))
            mock_executor.assert_called_with(
                *cmds, cpus=4, job_name='bwa_mem', mem=24, working_dir='tests/assets/jobs/testrun', log_commands=False
            )


class TestSamtoolsStatsMulti(TestPostDemultiplexing):
    def setUp(self):
        self.setup_stage(dm.SamtoolsStatsMulti)
        self.setup_post_alignment()

    def test_run(self):
        with patch_executor as mock_executor:
            assert self.stage._run() == 0
            cmds = []
            for run_element in self.run_elements:
                cmds.append(bash_commands.samtools_stats(
                    self.stage.bam_path(run_element),
                    self.stage.final_fastq_base(run_element) + '_samtools_stats.txt'
                ))
            mock_executor.assert_called_with(
                *cmds,
                cpus=1,
                job_name='samtoolsstats',
                mem=8,
                working_dir='tests/assets/jobs/testrun',
                log_commands=False
            )


class TestSamtoolsDepthMulti(TestPostDemultiplexing):
    def setUp(self):
        self.setup_stage(dm.SamtoolsDepthMulti)
        self.setup_post_alignment()

    def test_run(self):
        with patch_executor as mock_executor:
            assert self.stage._run() == 0
            cmds = []
            for run_element in self.run_elements:
                cmds.append(bash_commands.samtools_depth_command(
                    self.stage.job_dir,
                    self.stage.bam_path(run_element),
                    self.stage.final_fastq_base(run_element) + '_samtools.depth'
                ))
            mock_executor.assert_called_with(
                *cmds,
                cpus=1,
                job_name='samtoolsdepth',
                mem=6,
                working_dir='tests/assets/jobs/testrun',
                log_commands=False
            )


class TestPicardMarkDuplicateMulti(TestPostDemultiplexing):
    def setUp(self):
        self.setup_stage(dm.PicardMarkDuplicateMulti)
        self.setup_post_alignment()

    def test_run(self):
        with patch_executor as mock_executor:
            assert self.stage._run() == 0
            cmds = []
            for run_element in self.run_elements:
                cmds.append(bash_commands.picard_mark_dup_command(
                    self.stage.bam_path(run_element),
                    self.stage.bam_path(run_element)[:-len('.bam')] + '_markdup.bam',
                    self.stage.final_fastq_base(run_element) + '_markdup.metrics'
                ))
            mock_executor.assert_called_with(
                *cmds, cpus=1, job_name='picardMD', mem=12, working_dir='tests/assets/jobs/testrun'
            )


class TestPicardInsertSizeMulti(TestPostDemultiplexing):
    def setUp(self):
        self.setup_stage(dm.PicardInsertSizeMulti)
        self.setup_post_alignment()

    def test_run(self):
        with patch_executor as mock_executor:
            assert self.stage._run() == 0
            cmds = []
            for run_element in self.run_elements:
                cmds.append(bash_commands.picard_insert_size_command(
                    self.stage.bam_path(run_element),
                    self.stage.final_fastq_base(run_element) + '_insertsize.metrics',
                    self.stage.final_fastq_base(run_element) + '_insertsize.pdf'
                ))
            mock_executor.assert_called_with(
                *cmds, cpus=1, job_name='picardIS', mem=12, working_dir='tests/assets/jobs/testrun'
            )


class TestPicardGCBias(TestPostDemultiplexing):
    def setUp(self):
        self.setup_stage(dm.PicardGCBias)
        self.setup_post_alignment()

    @patch_executor
    def test_run(self, mocked_executor):
        self.stage.dataset.reference_genome.return_value = 'a_genome.fa'
        with patch('analysis_driver.pipelines.demultiplexing.bash_commands.picard_gc_bias', return_value='a_cmd'):
            self.stage._run()

        mocked_executor.assert_called_with(
            'a_cmd',
            'a_cmd',
            job_name='gc_bias',
            working_dir='tests/assets/jobs/testrun',
            cpus=1,
            mem=4
        )


class TestRunReview(TestAnalysisDriver):
    @patch('analysis_driver.pipelines.demultiplexing.rest_communication.post_entry')
    def test_run_review(self, mocked_post):
        r = dm.RunReview(dataset=NamedMock(real_name='test_dataset'))
        assert r._run() == 0
        mocked_post.assert_called_with(
            'actions', {'action_type': 'automatic_run_review', 'run_id': 'test_dataset'}, use_data=True
        )
