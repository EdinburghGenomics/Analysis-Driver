from unittest.mock import patch

from analysis_driver.reader.version_reader import Version_reader, get_versions
from analysis_driver.config import default as cfg
from tests.test_analysisdriver import TestAnalysisDriver


class TestVersion(TestAnalysisDriver):

    def test_get_version(self ):
        with patch('analysis_driver.reader.version_reader.Version_reader._get_stdout_from_command', return_value='1.2'):
            reader = Version_reader(tool_name='bwa', command=' -h | grep "Version')
            assert reader.get_version() == '1.2'

        reader = Version_reader(tool_name='bwa', command=None)
        assert reader.get_version() == '0.9'

    def test_get_stdout_from_command(self):
        reader = Version_reader(tool_name='echo', command=' 123')
        reader._get_stdout_from_command() == '123'


    @patch('analysis_driver.reader.version_reader.Version_reader._get_stdout_from_command', return_value='1.2')
    def test_get_versions(self, mocked_get_stdout):
        versions = get_versions()
        expected_keys = ['bamtools', 'bcbio', 'bcftools', 'bcl2fastq', 'bwa', 'fastqc',
                         'fastqscreen', 'gatk', 'sambamba', 'samblaster', 'samtools',
                         'seqtk', 'sickle', 'tabix', 'verifybamid', 'well_duplicate']
        assert sorted(versions.keys()) == sorted(expected_keys)


