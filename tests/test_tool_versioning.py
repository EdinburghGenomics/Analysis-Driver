from os.path import join
from unittest.mock import patch
from analysis_driver.tool_versioning import Toolset
from analysis_driver.config import Configuration
from tests import TestAnalysisDriver


def fake_check_version(self, toolname, executable, exp_version, version_cmd, dependency):
    return executable.split('_')[1] == str(exp_version)


class TestToolset(TestAnalysisDriver):
    def setUp(self):
        self.toolset = Toolset()
        self.toolset.versioning_cfg = Configuration(join(self.assets_path, 'tool_versioning.yaml'))
        self.toolset.select_type('non_human_sample_processing')

    @patch('analysis_driver.tool_versioning.Toolset.check_version', new=fake_check_version)
    def test_tool_paths(self):
        self.toolset.select_version(0)
        assert self.toolset['bwa'] == 'path/to/bwa_1.0'
        assert self.toolset['gatk'] == 'path/to/gatk_3'
        assert self.toolset['fastqc'] == 'path/to/fastqc_v0.11.5'
        assert self.toolset['samtools'] == 'path/to/samtools_1.3.1'
        assert self.toolset['qsub'] == '/bin/sh'
        assert 'qsub' in self.toolset.unversioned_tools
        assert self.toolset.tool_versions == {'bwa': 1.0, 'fastqc': 'v0.11.5', 'gatk': 3, 'samtools': '1.3.1', 'java': 8}

        self.toolset.select_version(1)
        assert self.toolset['bwa'] == 'path/to/bwa_1.1'
        assert self.toolset['gatk'] == 'path/to/gatk_4'
        assert self.toolset['fastqc'] == 'path/to/fastqc_v0.11.5'
        assert self.toolset['java'] == 'path/to/java_8'
        assert 'fastqc' in self.toolset.unversioned_tools
        assert self.toolset.tool_versions == {'bwa': 1.1, 'gatk': 4, 'samtools': '1.3.1', 'java': 8}

    @patch('analysis_driver.tool_versioning.Toolset.check_version')
    def test_tool_config_inheritance(self, mocked_check):
        self.toolset.select_version(1)
        mocked_check.assert_any_call(
            'bwa',
            'path/to/bwa_1.0',  # a bit wonky here since we've patched check_version, but the main test still stands
            1.1,
            'grep_version',
            None
        )
        mocked_check.assert_any_call(
            'gatk',
            'path/to/gatk_3',  # same here
            4,
            '{dependency} -jar {executable} -v 2>&1 | grep "GATK" | cut -d " " -f 2',
            'path/to/java_8'
        )
        mocked_check.assert_any_call(
            'samtools',
            'path/to/samtools_1.3.1',
            '1.3.1',
            'grep_version',
            None
        )
        mocked_check.assert_any_call(
            'java',
            'path/to/java_8',
            8,
            "{executable} -version 2>&1 | head -n 1 | cut -d ' ' -f 3 | sed 's/\"//g'",
            None
        )

    @patch('analysis_driver.tool_versioning.Toolset._get_stdout', return_value='1.0')
    def test_check_tool_version(self, mocked_stdout):
        assert self.toolset.check_version('bwa', 'path/to/bwa_1.0', 1.0, 'bwa -v')
        assert not self.toolset.check_version('bwa', 'path/to/bwa_1.0', 1.1, 'grep_version')
        mocked_stdout.assert_called_with(
            'path/to/bwa_1.0 2>&1 | grep "Version" | cut -d " " -f 2'
        )
