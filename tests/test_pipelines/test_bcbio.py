from unittest.mock import Mock, patch, mock_open

import yaml

from analysis_driver.pipelines.bcbio import FixUnmapped, BCBio
from analysis_driver.pipelines.demultiplexing import FastqFilter
from tests.test_analysisdriver import TestAnalysisDriver, NamedMock



class TestBCBio(TestAnalysisDriver):
    def test_run(self):
        dataset = NamedMock(
            real_name='test',
            user_sample_id='usertest',
            genome_version='hg38'
        )
        b = BCBio(dataset=dataset)

        patch_executor = patch('analysis_driver.pipelines.bcbio.executor.execute')
        patch_get_sample = patch('analysis_driver.pipelines.bcbio.clarity.get_sample', return_value=Mock(udf={
            'Analysis Type': 'gatk'
        }))
        patch_chdir = patch('os.chdir')
        yaml_content = {'fc_name': 'fc1'}
        # FIXME: have to patch the builtins to get 3.4 support patch directly in the file in 3.6
        # patch_open = patch('analysis_driver.pipelines.bcbio.open', new=mock_open(read_data=yaml.safe_dump(yaml_content)))
        patch_open = patch('builtins.open', new=mock_open(read_data=yaml.safe_dump(yaml_content)))

        with patch_executor as pexecute, patch_get_sample as pgetsample, patch_chdir, patch_open as popen:
            b._run()
            pgetsample.assert_called_once_with('test')
            bcb_cmd = (
                'path/to/bcbio/bin/bcbio_nextgen.py '
                'tests/assets/jobs/test/samples_test-merged/config/samples_test-merged.yaml '
                '-n 16 --workdir tests/assets/jobs/test/samples_test-merged/work'
            )
            # test first bcbio command's first argument which is just to command
            pexecute.mock_calls[0][0] == ''.join(bcb_cmd)
            bcb_cmd = (
                'path/to/bcbio/bin/bcbio_nextgen.py '
                'tests/assets/jobs/test/samples_test-merged/config/samples_test-merged.yaml '
                '-n 16 --workdir tests/assets/jobs/test/samples_test-merged/work'
            )
            # test last bcbio command's first argument which is just to command
            pexecute.mock_calls[3][0] == ''.join(bcb_cmd)
            # assert reading
            popen.assert_any_call('tests/assets/jobs/test/samples_test-merged/config/samples_test-merged.yaml', 'r')
            # assert writing
            popen.assert_called_with('tests/assets/jobs/test/samples_test-merged/config/samples_test-merged.yaml', 'w')
            popen().write.assert_called_with('fc_name: usertest\n')


class TestFixUnmapped(TestAnalysisDriver):
    def test_run(self):

        dataset = NamedMock(
            real_name='test',
            user_sample_id='usertest'
        )
        f = FixUnmapped(dataset=dataset)

        patch_executor = patch('analysis_driver.pipelines.bcbio.executor.execute')

        with patch_executor as pexecute:
            f._run()
            cmd = (
                'path/to/fix_dup_unmapped '
                 '-i tests/assets/jobs/test/samples_test-merged/final/usertest/usertest-ready.bam '
                 '-o tests/assets/jobs/test/samples_test-merged/final/usertest/usertest-ready_fixed.bam'
            )
            assert pexecute.call_args[0][0] == ''.join(cmd)
