import os
from shutil import rmtree
from unittest.mock import Mock, patch
from tests.test_analysisdriver import TestAnalysisDriver, NamedMock
from analysis_driver.config import output_file_config, cfg
from analysis_driver.pipelines import common


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
        for k in output_file_config['non_human_qc']:
            f = os.path.join(
                self.pseudo_links,
                output_file_config.output_dir_file('non_human_qc', k).format(
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
        patched_find_project = patch('egcg_core.clarity.find_project_name_from_sample', return_value='proj_' + self.sample_id)

        with patched_find_project, patched_archive, patched_crawler, patched_execute:
            stage = common.SampleDataOutput(dataset=self.dataset, output_fileset='non_human_qc')
            exit_status = stage.output_data(self.pseudo_links)

        output_files = os.path.join(cfg['sample']['output_dir'], 'proj_' + self.sample_id, self.sample_id)
        expected_outputs = ['samtools_stats.txt', 'taxa_identified.json', 'test_dataset.depth',
                            'test_dataset_R1_fastqc.html', 'test_dataset_R1_fastqc.zip',
                            'test_dataset_R1_screen.txt', 'test_dataset_R2_fastqc.html',
                            'test_dataset_R2_fastqc.zip', 'test_dataset_filter_snp.vcf.stats']

        o = sorted(os.listdir(output_files))
        assert exit_status == 0
        assert o == expected_outputs


class TestMergeFastqs(TestCommon):
    def setUp(self):
        super().setUp()
        self.stage = common.MergeFastqs(dataset=self.dataset)
        self.bcbio_csv = os.path.join(self.assets_path, 'samples_test_dataset.csv')
        self.job_dir = os.path.abspath(self.stage.job_dir)
        os.makedirs(self.job_dir, exist_ok=True)

    def tearDown(self):
        if os.path.isfile(self.bcbio_csv):
            os.remove(self.bcbio_csv)
        if os.path.isdir(self.job_dir):
            rmtree(self.job_dir)

    @patch('egcg_core.util.find_fastqs')
    def test_find_fastqs(self, mocked_find):
        fake_run_element = {'run_id': 'a_run', 'project_id': 'a_project', 'lane': 1}
        self.stage._find_fastqs_for_run_element(fake_run_element)
        mocked_find.assert_called_with(os.path.join(cfg['sample']['input_dir'], 'a_run'), 'a_project', 'test_dataset', 1)

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

    def test_run_overwrite(self):
        infastqs = ['fastq1_R1', 'fastq1_R2', 'fastq2_R1', 'fastq2_R2']
        with patch.object(common.MergeFastqs, 'find_fastqs_for_sample', return_value=infastqs), \
                patch('egcg_core.executor.execute') as mock_execute:
            mock_execute().join.return_value = 0
            merge_dir = os.path.join(self.job_dir, 'merged')
            os.makedirs(merge_dir)
            self._touch(os.path.join(merge_dir, 'user_sample_id_R1.fastq.gz'))
            self._touch(os.path.join(merge_dir, 'user_sample_id_R2.fastq.gz'))
            self.dataset.user_sample_id = 'user_sample_id'

            assert os.path.exists(os.path.join(merge_dir, 'user_sample_id_R1.fastq.gz'))
            assert os.path.exists(os.path.join(merge_dir, 'user_sample_id_R2.fastq.gz'))

            self.stage._run()

            # overwrite has deleted the files
            assert not os.path.exists(os.path.join(merge_dir, 'user_sample_id_R1.fastq.gz'))
            assert not os.path.exists(os.path.join(merge_dir, 'user_sample_id_R2.fastq.gz'))
            command = 'path/to/bcbio/bin/bcbio_prepare_samples.py ' \
                      '--out tests/assets/jobs/test_dataset/merged' \
                      ' --csv tests/assets/jobs/test_dataset/samples_test_dataset.csv'
            mock_execute.assert_any_call(command, job_name='bcbio_prepare_samples', working_dir='tests/assets/jobs/test_dataset'),


class TestSampleReview(TestCommon):
    @patch('analysis_driver.pipelines.common.rest_communication.post_entry')
    def test_sample_review(self, mocked_post):
        r = common.SampleReview(dataset=self.dataset)
        assert r._run() == 0
        mocked_post.assert_called_with(
            'actions', {'action_type': 'automatic_sample_review', 'sample_id': 'test_dataset'}, use_data=True
        )
