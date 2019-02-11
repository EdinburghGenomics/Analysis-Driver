import os
from shutil import rmtree
from unittest.mock import patch
from tests.test_analysisdriver import TestAnalysisDriver, NamedMock
from analysis_driver.pipelines import projects

test_projects = os.path.join(TestAnalysisDriver.assets_path, 'test_projects')
relatedness_outfiles = os.path.join(test_projects, 'relatedness_outfiles')


class TestProjects(TestAnalysisDriver):
    @staticmethod
    def fake_dataset(sample_ids, pipeline_used=None):
        samples_processed = []
        for s in sample_ids:
            sample = {'sample_id': s, 'user_sample_id': 'uid_' + s}
            if pipeline_used:
                sample['aggregated'] = {'most_recent_proc': {'pipeline_used': {'name': pipeline_used}}}
            samples_processed.append(sample)

        return NamedMock(
            real_name='test_dataset',
            samples_processed=samples_processed
        )

    def test_2_samples(self):
        pipeline = projects.build_pipeline(self.fake_dataset(['10015AT0004', '10015AT0003'], 'bcbio'))
        assert len(pipeline.previous_stages) > 0

    def test_1_sample(self):
        pipeline = projects.build_pipeline(self.fake_dataset(['10015AT0004']))
        assert len(pipeline.previous_stages) == 0

    def test_non_human(self):
        pipeline = projects.build_pipeline(self.fake_dataset(['sample1', 'sample2', 'sample3'], 'qc'))
        assert len(pipeline.previous_stages) == 0

    def test_combinegvcf_batching(self):
        sample_ids = ['sample' + str(n) for n in range(1, 31)]
        d = self.fake_dataset(sample_ids, 'bcbio')

        with patch('analysis_driver.pipelines.projects.find_file', return_value='a_file'):
            final_stage = projects.build_pipeline(d)

        parse = final_stage.previous_stages[0].previous_stages[0].previous_stages[0]
        genotype_gvcfs = parse.previous_stages[0].previous_stages[0]
        assert genotype_gvcfs.previous_stages[0].gvcfs == tuple('a_file' for _ in range(25))  # batch 1
        assert genotype_gvcfs.previous_stages[1].gvcfs == tuple('a_file' for _ in range(5))  # batch 2


class TestMD5Sum(TestAnalysisDriver):
    def setUp(self):
        os.makedirs(relatedness_outfiles, exist_ok=True)
        self.filename = os.path.join(relatedness_outfiles, 'an_output_file.txt')
        open(self.filename, 'w').close()

    def tearDown(self):
        os.remove(self.filename)
        rmtree(relatedness_outfiles)

    @patch.object(projects.MD5Sum, 'job_dir', new=test_projects)
    @patch('egcg_core.executor.execute')
    def test_run(self, mocked_execute):
        md5 = projects.MD5Sum(dataset=NamedMock(real_name='test_dataset'))
        with patch('os.makedirs'):
            md5._run()

        mocked_execute.assert_called_with(
            'path/to/md5sum %s > %s' % (self.filename, self.filename + '.md5'),
            cpus=1,
            mem=2,
            job_name='md5sum',
            log_commands=False,
            working_dir=test_projects
        )


class TestOutput(TestAnalysisDriver):
    def setUp(self):
        self.o = projects.Output(dataset=NamedMock(real_name='test_dataset'))

    @patch.object(projects.toolset, 'write_to_yaml')
    @patch.object(projects, 'create_output_links')
    @patch.object(projects, 'output_data_and_archive')
    @patch.object(projects, 'output_file_config')
    def test_run(self, mocked_outfile_config, mocked_output_archive, mocked_output_links, mocked_write):
        with patch('analysis_driver.segmentation.BasicStage.job_dir', new=test_projects):
            self.o._run()
            mocked_output_archive.assert_called_with(relatedness_outfiles, 'path/to/input/dir/test_dataset')
            mocked_outfile_config.set_pipeline_type.assert_called_with('project_process')
            mocked_write.assert_called_with(os.path.join(relatedness_outfiles, 'program_versions.yaml'))
            mocked_output_links.assert_called_with(
                test_projects,
                mocked_outfile_config,
                relatedness_outfiles,
                project_id='test_dataset'
            )
