import os
import multiprocessing
from unittest.mock import patch
from analysis_driver.tool_versioning import Toolset
from analysis_driver.config import Configuration
from tests import TestAnalysisDriver


def fake_check_version(self, toolname, executable, exp_version, version_cmd, dependency):
    return executable.split('_')[1] == str(exp_version)


class TestToolset(TestAnalysisDriver):
    def setUp(self):
        self.toolset = Toolset()
        self.toolset.versioning_cfg = Configuration(os.path.join(self.assets_path, 'tool_versioning.yaml'))
        self.program_versions = os.path.join(self.assets_path, 'program_versions.yaml')
        self.toolset.configure('non_human_sample_processing', 0, self.program_versions)

    def tearDown(self):
        if os.path.isfile(self.program_versions):
            os.remove(self.program_versions)

    @patch('analysis_driver.tool_versioning.Toolset.check_version', new=fake_check_version)
    def test_tool_paths(self):
        assert self.toolset['bwa'] == 'path/to/bwa_1.0'
        assert self.toolset['gatk'] == 'path/to/gatk_3'
        assert self.toolset['fastqc'] == 'path/to/fastqc_v0.11.5'
        assert self.toolset['samtools'] == 'path/to/samtools_1.3.1'
        assert self.toolset['qsub'] == '/bin/sh'
        assert self.toolset.tool_versions == {'bwa': 1.0, 'fastqc': 'v0.11.5', 'gatk': 3, 'samtools': '1.3.1', 'java': 8}

        self.toolset.configure('non_human_sample_processing', 1, self.program_versions)
        assert self.toolset['bwa'] == 'path/to/bwa_1.1'
        assert self.toolset['gatk'] == 'path/to/gatk_4'
        assert self.toolset['fastqc'] == 'path/to/fastqc_v0.11.5'
        assert self.toolset['samtools'] == 'path/to/samtools_1.3.1'
        assert self.toolset['java'] == 'path/to/java_8'
        assert 'fastqc' in self.toolset.unversioned_tools  # no longer versioned
        assert self.toolset.tool_versions == {
            'bwa': 1.1,  # updated
            'gatk': 4,  # updated, new version_cmd
            'samtools': '1.3.1',  # inherited
            'java': 8  # inherited
        }

    @patch('analysis_driver.tool_versioning.Toolset.check_version')
    def test_tool_config_inheritance(self, mocked_check):
        self.toolset.configure('non_human_sample_processing', 1, self.program_versions)
        for t in ('bwa', 'gatk', 'samtools', 'java'):
            _ = self.toolset[t]

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

    @patch('analysis_driver.tool_versioning.Toolset.check_version', new=fake_check_version)
    def test_update_versions_file(self):
        self.toolset.update_versions_file('fastqc', 'v0.11.5')
        with open(self.program_versions, 'r') as f:
            assert f.read() == 'fastqc: v0.11.5\n'

        self.toolset.update_versions_file('gatk', 3)
        with open(self.program_versions, 'r') as f:
            assert f.read() == 'fastqc: v0.11.5\ngatk: 3\n'

    @patch('analysis_driver.tool_versioning.Toolset.check_version', new=fake_check_version)
    def test_multiprocessing(self):
        def run_toolset():
            _ = self.toolset['fastqc']  # reference a tool - should appear in versions_file
            _ = self.toolset['bcbio']  # unversioned, should not appear below
            _ = self.toolset['gatk']  # should appear below along with Java

        procs = [multiprocessing.Process(target=run_toolset) for _ in range(3)]
        for p in procs:
            p.start()

        for p in procs:
            p.join()

        with open(self.program_versions, 'r') as f:
            assert f.read() == 'fastqc: v0.11.5\ngatk: 3\njava: 8\n'
