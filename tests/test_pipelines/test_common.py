import os
from shutil import rmtree
from unittest.mock import Mock, patch
from tests.test_analysisdriver import TestAnalysisDriver, NamedMock
from analysis_driver.config import OutputFileConfiguration
from analysis_driver.pipelines import common
from analysis_driver.config import cfg


class TestCommon(TestAnalysisDriver):
    def setUp(self):
        self.sample_id = 'test_dataset'
        self.dataset = NamedMock(real_name=self.sample_id)


class TestSampleDataOutput(TestCommon):
    data_output = os.path.join(TestAnalysisDriver.assets_path, 'data_output')

    def setUp(self):
        super().setUp()
        self.pseudo_links = os.path.join(self.data_output, 'pseudo_links')
        self.to_dir = os.path.join(self.data_output, 'to', '')
        os.makedirs(self.pseudo_links, exist_ok=True)
        output_cfg = OutputFileConfiguration('non_human_qc')
        for k in output_cfg.content:
            f = os.path.join(
                self.pseudo_links,
                output_cfg.output_dir_file(k).format(
                    sample_id=self.sample_id, user_sample_id=self.sample_id
                )
            )
            open(f, 'a').close()

    def tearDown(self):
        for d in (self.pseudo_links, self.to_dir):
            if os.path.isdir(d):
                rmtree(d)

    def test_output_data(self):
        patched_archive = patch('analysis_driver.transfer_data.archive_management.archive_directory', return_value=True)
        patched_crawler = patch('analysis_driver.pipelines.common.SampleCrawler')
        patched_execute = patch('egcg_core.executor.execute', return_value=Mock(join=Mock(return_value=0)))
        patched_cfg = patch('analysis_driver.pipelines.common.cfg', new={'output_dir': self.to_dir})
        patched_find_project = patch('egcg_core.clarity.find_project_name_from_sample', return_value='proj_' + self.sample_id)

        with patched_find_project, patched_archive, patched_crawler, patched_execute, patched_cfg:
            stage = common.SampleDataOutput(dataset=self.dataset, output_fileset='non_human_qc')
            exit_status = stage.output_data(self.pseudo_links)

        output_files = os.path.join(self.to_dir, 'proj_' + self.sample_id, self.sample_id)
        expected_outputs = ['samtools_stats.txt', 'taxa_identified.json', 'test_dataset.depth',
                            'test_dataset_R1_fastqc.html', 'test_dataset_R1_fastqc.zip',
                            'test_dataset_R1_screen.txt', 'test_dataset_R2_fastqc.html',
                            'test_dataset_R2_fastqc.zip']

        o = sorted(os.listdir(output_files))
        assert exit_status == 0
        assert o == expected_outputs


class TestMergeFastqs(TestCommon):
    def setUp(self):
        super().setUp()
        self.stage = common.MergeFastqs(dataset=self.dataset)
        self.bcbio_csv = os.path.join(self.assets_path, 'samples_test_dataset.csv')

    def tearDown(self):
        if os.path.isfile(self.bcbio_csv):
            os.remove(self.bcbio_csv)

    @patch('egcg_core.util.find_fastqs')
    def test_find_fastqs(self, mocked_find):
        fake_run_element = {'run_id': 'a_run', 'project_id': 'a_project', 'lane': 1}
        self.stage._find_fastqs_for_run_element(fake_run_element)
        mocked_find.assert_called_with(os.path.join(cfg['input_dir'], 'a_run'), 'a_project', 'test_dataset', 1)

    def test_write_bcbio_csv(self):
        self.dataset.user_sample_id = 'a_user_sample_id'
        patched_job_dir = patch('analysis_driver.pipelines.common.MergeFastqs.job_dir', new=self.assets_path)
        with patched_job_dir:
            self.stage._write_bcbio_csv(['test_R1.fastq', 'test_R2.fastq'])

        with open(self.bcbio_csv) as f:
            content = f.read()
            assert content == (
                'samplename,description\n'
                'test_R1.fastq,a_user_sample_id\n'
                'test_R2.fastq,a_user_sample_id\n'
            )
