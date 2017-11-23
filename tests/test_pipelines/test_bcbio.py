import yaml
from unittest.mock import Mock, patch, mock_open
from tests.test_analysisdriver import TestAnalysisDriver, NamedMock
from analysis_driver.pipelines.bcbio import FixUnmapped, BCBio
from analysis_driver.config import etc_config


class TestBCBio(TestAnalysisDriver):
    def test_run(self):
        dataset = NamedMock(
            real_name='test',
            user_sample_id='usertest',
            genome_version='hg38'
        )
        b = BCBio(dataset=dataset)

        patch_executor = patch('analysis_driver.pipelines.bcbio.executor.execute')
        patch_get_sample = patch(
            'analysis_driver.pipelines.bcbio.clarity.get_sample', return_value=Mock(udf={'Analysis Type': 'gatk'})
        )
        patch_chdir = patch('os.chdir')
        yaml_content = {'fc_name': 'fc1'}
        # FIXME: have to patch the builtins to get 3.4 support patch directly in the file in 3.6
        patch_open = patch('builtins.open', new=mock_open(read_data=yaml.safe_dump(yaml_content)))

        with patch_executor as pexecute, patch_get_sample as pgetsample, patch_chdir, patch_open as popen:
            b._run()
            pgetsample.assert_called_once_with('test')

            obs = pexecute.call_args_list[0][0][0]
            exp = ' '.join(
                (
                    'path/to/bcbio/bin/bcbio_nextgen.py -w template',
                    etc_config('bcbio_alignment_hg38_gatk.yaml'),
                    'path/to/jobs/test/samples_test-merged path/to/jobs/test/samples_test-merged.csv'
                )
            )
            assert obs == exp

            assert pexecute.call_args_list[1][0][0] == (
                'path/to/bcbio/bin/bcbio_nextgen.py '
                'path/to/jobs/test/samples_test-merged/config/samples_test-merged.yaml '
                '-n 16 --workdir path/to/jobs/test/samples_test-merged/work'
            )

            # assert reading
            popen.assert_any_call('path/to/jobs/test/samples_test-merged/config/samples_test-merged.yaml', 'r')
            # assert writing
            popen.assert_called_with('path/to/jobs/test/samples_test-merged/config/samples_test-merged.yaml', 'w')
            popen().write.assert_called_with('fc_name: usertest\n')


class TestFixUnmapped(TestAnalysisDriver):
    @patch('analysis_driver.pipelines.bcbio.executor.execute')
    def test_run(self, mocked_execute):
        f = FixUnmapped(
            dataset=NamedMock(
                real_name='test',
                user_sample_id='usertest'
            )
        )
        f._run()
        cmd = (
            'path/to/fix_dup_unmapped '
            '-i path/to/jobs/test/samples_test-merged/final/usertest/usertest-ready.bam '
            '-o path/to/jobs/test/samples_test-merged/final/usertest/usertest-ready_fixed.bam'
        )
        assert mocked_execute.call_args[0][0] == ''.join(cmd)
