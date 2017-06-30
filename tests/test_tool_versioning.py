from unittest.mock import patch
from analysis_driver.tool_versioning import toolset
from analysis_driver.config import cfg, _etc_config
cfg.load_config_file(_etc_config('example_analysisdriver.yaml'))


def fake_check_version(toolname, executable, exp_version, version_cmd):
    return executable.split('_')[1] == str(exp_version)


def test_tool_versioning():
    toolset.check_version = fake_check_version

    toolset.select_toolset(0)
    assert toolset['bwa'] == 'path/to/bwa_1.0'
    assert toolset['gatk'] == 'path/to/GenomeAnalysisTK.jar_3'
    assert toolset['bcl2fastq'] == 'path/to/bcl2fastq_1.0.4'
    assert toolset['qsub'] == '/bin/sh'

    toolset.select_toolset(1)
    assert toolset['bwa'] == 'path/to/bwa_1.1'
    assert toolset['gatk'] == 'path/to/GenomeAnalysisTK.jar_4'
    assert toolset['bcl2fastq'] == 'path/to/bcl2fastq_1.0.4'
    assert toolset['qsub'] == '/bin/sh'


@patch('analysis_driver.tool_versioning.Toolset._get_stdout')
def test_check_tool_version(mocked_stdout):
    mocked_stdout.return_value = '1.0'
    assert toolset.check_version('bwa', 'path/to/bwa_1.0', 1.0, 'bwa -v')
    assert not toolset.check_version('bwa', 'path/to/bwa_1.0', 1.1, 'bwa -v')
