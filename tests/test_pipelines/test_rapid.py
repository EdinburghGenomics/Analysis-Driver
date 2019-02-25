import os
from shutil import rmtree
from egcg_core.config import cfg
from unittest.mock import Mock, patch
from analysis_driver.pipelines import rapid
from tests.test_analysisdriver import TestAnalysisDriver, NamedMock


test_run_dir = os.path.join(TestAnalysisDriver.assets_path, 'rapid_analysis', 'a_run')
patched_uid = patch('egcg_core.clarity.get_user_sample_name', return_value='a_uid')
patched_executor = patch(
    'analysis_driver.pipelines.demultiplexing.executor.execute',
    return_value=Mock(join=Mock(return_value=0))
)


class TestDragen(TestAnalysisDriver):
    @patched_uid
    @patched_executor
    def test_run(self, mocked_execute, mocked_uid):
        d = rapid.Dragen(
            dataset=NamedMock(
                real_name='a_run',
                type='run',
                rapid_samples_by_lane={'2': NamedMock(real_name='a_sample')},
                sample_sheet_file='path/to/SampleSheet_analysis_driver.csv'
            )
        )
        d._run()
        mocked_execute.assert_called_with(
            'if [ -d /path/to/dragen_staging/a_run_2 ]; then rm -r /path/to/dragen_staging/a_run_2; fi; '
            'mkdir -p /path/to/dragen_staging/a_run_2; dragen -r /path/to/dragen_reference '
            '--bcl-input-dir tests/assets/data_transfer/from/a_run --bcl-only-lane 2 '
            '--output-directory /path/to/dragen_staging/a_run_2 --output-file-prefix a_uid '
            '--enable-variant-caller true --vc-sample-name a_sample '
            '--bcl-sample-sheet path/to/SampleSheet_analysis_driver.csv; '
            'rsync -rLD /path/to/dragen_staging/a_run_2/ tests/assets/jobs/a_run/rapid_analysis_2; '
            'rm -r /path/to/dragen_staging/a_run_2',
            prelim_cmds=['ulimit -n 65535', 'ulimit -c 0', 'ulimit -s 8192', 'ulimit -u 16384', 'ulimit -i 1029522'],
            job_name='dragen',
            queue='dragen',
            working_dir='tests/assets/jobs/a_run',
            exclusive=True
        )


class TestDragenMetrics(TestAnalysisDriver):
    @patch('analysis_driver.pipelines.rapid.rest_communication.post_or_patch')
    def test_run(self, mocked_post_or_patch):
        d = rapid.DragenMetrics(
            dataset=NamedMock(
                real_name='a_run',
                type='run',
                rapid_samples_by_lane={'2': NamedMock(real_name='a_sample'), '4': NamedMock(real_name='another_sample')}
            )
        )

        with patch.object(d.__class__, 'job_dir', new=test_run_dir):
            d._run()

        mocked_post_or_patch.assert_called_with(
            'samples',
            [
                {
                    'sample_id': 'a_sample',
                    'rapid_metrics': {
                        'var_calling': {'ti_tv_ratio': 1.9, 'het_hom_ratio': 1.53},
                        'mapping': {
                            'total_reads': 911979369,
                            'duplicate_reads': 107673375,
                            'pc_duplicates': 11.76,
                            'mapped_reads': 872974134,
                            'pc_mapped': 95.67,
                            'unique_mapped_reads': 765300754,
                            'pc_unique_mapped': 83.87,
                            'mean_coverage': 35.52,
                            'median_coverage': 35.9,
                            'genome_size': 3217346912,
                            'pc_genome_with_coverage': {
                                '0-1': 7.43,
                                '1-2': 0.31,
                                '2-5': 0.5,
                                '5-10': 0.61,
                                '10-20': 4.21,
                                '20-30': 12.82,
                                '30-40': 38.61,
                                '40-50': 29.09,
                                '50-100': 5.85,
                                '100+': 0.08
                            }
                        }
                    }
                },
                {
                    'sample_id': 'another_sample',
                    'rapid_metrics': {
                        'var_calling': {'ti_tv_ratio': 1.94, 'het_hom_ratio': 1.54},
                        'mapping': {
                            'total_reads': 911979370,
                            'duplicate_reads': 107673376,
                            'pc_duplicates': 11.77,
                            'mapped_reads': 872974135,
                            'pc_mapped': 95.68,
                            'unique_mapped_reads': 765300755,
                            'pc_unique_mapped': 83.88,
                            'mean_coverage': 35.55,
                            'median_coverage': 35.94,
                            'genome_size': 3217346913,
                            'pc_genome_with_coverage': {
                                '0-1': 7.44,
                                '1-2': 0.32,
                                '2-5': 0.51,
                                '5-10': 0.62,
                                '10-20': 4.22,
                                '20-30': 12.83,
                                '30-40': 38.62,
                                '40-50': 29.1,
                                '50-100': 5.86,
                                '100+': 0.39
                            }
                        }
                    }
                }
            ],
            id_field='sample_id'
        )


class TestDragenOutput(TestAnalysisDriver):
    def setUp(self):
        self.dragen_output_dir = os.path.join(test_run_dir, 'rapid_analysis_2')
        self.staging_dir = os.path.join(test_run_dir, 'linked_output_files_2')
        self.sample_output_dir = os.path.join(cfg['sample']['output_dir'], 'a_rapid_project', 'a_sample')
        os.makedirs(self.staging_dir, exist_ok=True)
        os.makedirs(self.dragen_output_dir, exist_ok=True)

        self.dragen_output_files = [
            os.path.join(self.dragen_output_dir, 'a_uid.' + ext)
            for ext in ('vcf.gz', 'vcf.gz.tbi', 'vcf.gz.md5sum')
        ]

        for f in self.dragen_output_files:
            open(f, 'w').close()

    @patched_executor
    @patch('analysis_driver.pipelines.rapid.bash_commands.md5sum', return_value='an_md5_command')
    @patch('analysis_driver.transfer_data.archive_management.archive_directory', return_value=True)
    def test_output_data(self, mocked_archive, mocked_md5, mocked_executor):
        d = rapid.DragenOutput(
            dataset=NamedMock(
                real_name='a_run',
                type='run',
                rapid_samples_by_lane={'2': NamedMock(real_name='a_sample'), '4': NamedMock(real_name='another_sample')}
            )
        )

        for f in self.dragen_output_files:
            assert os.path.isfile(f)

        with patched_uid, patch.object(d.__class__, 'job_dir', new=test_run_dir):
            d.output_data(2, NamedMock(real_name='a_sample', project=NamedMock(real_name='a_rapid_project')))

        for f in self.dragen_output_files:
            assert not os.path.isfile(f)

        for ext in ('vcf.gz', 'vcf.gz.tbi', 'vcf.gz.md5'):
            assert os.path.isfile(os.path.join(self.sample_output_dir, 'a_uid.' + ext))

        assert mocked_md5.call_count == 2
        mocked_archive.assert_called_with(self.sample_output_dir)
        mocked_executor.assert_called_with(
            'an_md5_command',
            'an_md5_command',
            job_name='md5sum',
            working_dir=test_run_dir,
            cpus=1,
            mem=2,
            log_commands=False
        )

    def tearDown(self):
        for d in (self.sample_output_dir, self.staging_dir):
            if os.path.isdir(d):
                rmtree(d)
