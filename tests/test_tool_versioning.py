from unittest.mock import patch
from analysis_driver.tool_versioning import toolset
from analysis_driver.config import cfg, _etc_config
cfg.load_config_file(_etc_config('example_analysisdriver.yaml'))


@patch('analysis_driver.tool_versioning.Toolset.check_version')
def test_tool_versioning(mocked_check_version):
    toolset.select_toolset_version(0)
    assert toolset.tools == {
        'bwa': 'path/to/bwa',
        'gatk': 'path/to/GenomeAnalysisTK3.jar',
        'bcl2fastq': 'path/to/bcl2fastq'
    }

    toolset.select_toolset_version(1)
    assert toolset.tools == {
        'bwa': 'path/to/bwa',
        'gatk': 'path/to/GenomeAnalysisTK4.jar',
        'bcl2fastq': 'path/to/bcl2fastq'
    }


@patch('analysis_driver.tool_versioning.Toolset._get_stdout')
def test_check_tool_versio(mocked_stdout):
    mocked_stdout.return_value = '1.0'
    toolset.check_version('bwa', 1.0, 'bwa -v')
    assert toolset.toolset_valid

    toolset.check_version('bwa', 1.1, 'bwa -v')
    assert not toolset.toolset_valid
