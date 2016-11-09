from unittest.mock import patch
from analysis_driver.reader.version_reader import VersionReader, get_versions
from tests.test_analysisdriver import TestAnalysisDriver


class TestVersion(TestAnalysisDriver):
    def test_get_version(self):
        with patch('analysis_driver.reader.version_reader.VersionReader._get_stdout_from_command',
                   return_value='1.2'):
            reader = VersionReader(tool_name='bwa', command=' -h | grep "Version')
            assert reader.get_version() == '1.2'

        reader = VersionReader(tool_name='bwa', command=None)
        assert reader.get_version() == '0.9'

    def test_get_stdout_from_command(self):
        reader = VersionReader(tool_name='echo', command='{executable} 123')
        assert reader._get_stdout_from_command() == '123'

    def test_get_versions(self):
        with patch('analysis_driver.reader.version_reader.VersionReader.get_version', return_value='1.2'):
            observed = get_versions()

        expected_keys = ['bamtools', 'bcbio', 'bcftools', 'bcl2fastq', 'bwa', 'fastqc',
                         'fastqscreen', 'gatk', 'sambamba', 'samblaster', 'samtools',
                         'seqtk', 'tabix', 'verifybamid', 'well_duplicate', 'biobambam_sortmapdup']
        expected = dict((k, '1.2') for k in expected_keys)
        assert observed == expected
