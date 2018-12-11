import os
import shutil
from os.path import join
from unittest.mock import Mock, patch
from egcg_core.constants import ELEMENT_PROJECT_ID, ELEMENT_LANE, ELEMENT_SAMPLE_INTERNAL_ID
from tests.test_analysisdriver import TestAnalysisDriver, NamedMock
from analysis_driver.util import bash_commands
from analysis_driver.pipelines import demultiplexing as dm

patch_executor = patch(
    'analysis_driver.pipelines.demultiplexing.executor.execute',
    return_value=Mock(join=Mock(return_value=0))
)


class TestPhixDetection(TestAnalysisDriver):
    def test_run(self):
        dataset = NamedMock(
            real_name='test'
        )
        f = dm.PhixDetection(dataset=dataset)
        fake_fastq_pairs = [
            [('L1_R1_001.fastq.gz', 'L1_R2_001.fastq.gz')], [('L2_R1_001.fastq.gz', 'L2_R2_001.fastq.gz')]
        ]
        patch_run_crawler = patch('analysis_driver.pipelines.demultiplexing.RunCrawler', autospec=True)
        patch_find = patch('egcg_core.util.find_all_fastq_pairs', side_effect=fake_fastq_pairs)

        with patch_find, patch_executor as pexecute, patch_run_crawler as prun_crawler:
            f._run()

            assert pexecute.call_args[0][0] == (
                'set -o pipefail; path/to/bwa mem -t 16 /path/to/phix.fa L1_R1_001.fastq.gz | '
                'path/to/samtools view -F 4 | cut -f 1 | sort -u > L1_phix_read_name.list'
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
            f._run()

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

    def test_run(self):

        # Run info state 150 cycle for first read and 8 index cycle + 50 cycle for second read = 208
        run_info = Mock(reads=Mock(upstream_read=Mock(attrib={'NumCycles': '150'}), index_lengths=8))
        dataset = NamedMock(real_name='testrun', run_info=run_info)
        # cycle extracted states 310 cycle done
        pcycles = patch('analysis_driver.quality_control.interop_metrics.get_cycles_extracted',
                        return_value=range(1, 311))

        with pcycles as mcycle:
            # No sleep time
            self.stage = dm.WaitForRead2(dataset=dataset)
            self.stage._run()
            assert mcycle.call_count == 1

        pcycles = patch('analysis_driver.quality_control.interop_metrics.get_cycles_extracted',
                        side_effect=[range(1, 208), range(1, 209)])
        with pcycles as mcycle, patch('time.sleep') as msleep:
            self.stage = dm.WaitForRead2(dataset=dataset)
            self.stage._run()
            assert mcycle.call_count == 2
            assert msleep.call_count == 1


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
        os.makedirs(self.stage.fastq_intermetiate_dir, exist_ok=True)

        for run_element in self.run_elements:
            s_dir = join(self.stage.fastq_intermetiate_dir, run_element[ELEMENT_PROJECT_ID], run_element[ELEMENT_SAMPLE_INTERNAL_ID])
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
        patch_get_sample = patch('egcg_core.clarity.get_species_from_sample', return_value='Homo sapiens')
        with patch_get_sample, patch_executor as mock_executor:
            self.stage._run()
            cmds = []
            for run_element in self.run_elements:
                cmds.append(bash_commands.bwa_mem_biobambam(
                    self.stage.fastq_pair(run_element),
                    '/path/to/genome.fa',
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
            self.stage._run()
            cmds = []
            for run_element in self.run_elements:
                cmds.append(bash_commands.samtools_stats(
                    self.stage.bam_path(run_element),
                    self.stage.real_fastq_base(run_element) + '_samtools_stats.txt'
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
            self.stage._run()
            cmds = []
            for run_element in self.run_elements:
                cmds.append(bash_commands.samtools_depth_command(
                    self.stage.job_dir,
                    self.stage.bam_path(run_element),
                    self.stage.real_fastq_base(run_element) + '_samtools.depth'
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
            self.stage._run()
            cmds = []
            for run_element in self.run_elements:
                cmds.append(bash_commands.picard_mark_dup_command(
                    self.stage.bam_path(run_element),
                    self.stage.bam_path(run_element)[:-len('.bam')] + '_markdup.bam',
                    self.stage.real_fastq_base(run_element) + '_markdup.metrics'
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
            self.stage._run()
            cmds = []
            for run_element in self.run_elements:
                cmds.append(bash_commands.picard_insert_size_command(
                    self.stage.bam_path(run_element),
                    self.stage.real_fastq_base(run_element) + '_insertsize.metrics',
                    self.stage.real_fastq_base(run_element) + '_insertsize.pdf'
                ))
            mock_executor.assert_called_with(
                *cmds, cpus=1, job_name='picardIS', mem=12, working_dir='tests/assets/jobs/testrun'
            )


class TestRunReview(TestAnalysisDriver):
    @patch('analysis_driver.pipelines.demultiplexing.rest_communication.post_entry')
    def test_run_review(self, mocked_post):
        r = dm.RunReview(dataset=NamedMock(real_name='test_dataset'))
        assert r._run() == 0
        mocked_post.assert_called_with(
            'actions', {'action_type': 'automatic_run_review', 'run_id': 'test_dataset'}, use_data=True
        )
