import os
from unittest.mock import Mock, patch
from tests.test_analysisdriver import TestAnalysisDriver, NamedMock
from analysis_driver.pipelines import projects

test_projects = os.path.join(TestAnalysisDriver.assets_path, 'test_projects', 'test_dataset')
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


class TestOutput(TestAnalysisDriver):
    def setUp(self):
        self.o = projects.Output(dataset=NamedMock(real_name='test_dataset'))
        self.output_files = [
            os.path.join(test_projects, 'test_dataset' + e)
            for e in ('_genotype_gvcfs.vcf.gz', '.relatedness2', '.ped_check.csv', '.relatedness_output.egc', '.relatedness_output.gel')
        ]
        for f in self.output_files:
            open(f, 'w').close()

    def tearDown(self):
        for f in self.output_files:
            os.remove(f)

    @patch('egcg_core.executor.execute', return_value=Mock(join=Mock(return_value=0)))
    @patch.object(projects.toolset, 'write_to_yaml')
    @patch.object(projects, 'output_data_and_archive')
    def test_run(self, mocked_output_archive, mocked_write, mocked_execute):
        with patch('analysis_driver.segmentation.BasicStage.job_dir', new=test_projects):
            self.o._run()
            mocked_execute.assert_called_with(
                *[
                    'path/to/md5sum {f} > {f}.md5'.format(f=os.path.join(relatedness_outfiles, os.path.basename(f)))
                    for f in sorted(self.output_files)
                ],
                cpus=1,
                mem=2,
                job_name='md5sum',
                log_commands=False,
                working_dir=test_projects
            )
            mocked_output_archive.assert_called_with(relatedness_outfiles, 'path/to/input/dir/test_dataset')
            mocked_write.assert_called_with(os.path.join(relatedness_outfiles, 'program_versions.yaml'))
