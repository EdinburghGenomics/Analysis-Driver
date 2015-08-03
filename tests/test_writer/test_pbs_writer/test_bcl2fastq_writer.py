__author__ = 'mwham'
from tests.test_writer.test_pbs_writer.test_pbs_writer import TestPBSWriter
from analysis_driver.writer.pbs_writer.bcl2fastq_writer import BCL2FastqWriter
import os.path


class TestBCL2FastqWriter(TestPBSWriter):
    def setUp(self):
        super().setUp()
        self.writer = BCL2FastqWriter(
            self.tmp_script,
            walltime='24',
            cpus='2',
            mem='1',
            job_name='test',
            log_file='test.log',
            queue='test'
        )

    def test_bcl2fastq(self):
        self.writer._bcl2fastq(mask='this,that,other', input_dir=self.assets_path, fastq_path='a_fastq_path')
        expected_ending = 'bcl2fastq -l INFO --runfolder-dir %s --output-dir %s --sample-sheet %s ' \
                          '--use-bases-mask %s\n\n' % (self.assets_path,
                                                       'a_fastq_path',
                                                       os.path.join(self.assets_path, 'SampleSheet.csv'),
                                                       'this,that,other')
        print(self.writer.script)
        print(expected_ending)
        assert self.writer.script.endswith(expected_ending)

    def test_write(self):
        self.writer.write(
            job_dir=self.assets_path,
            mask='this,that,other',
            input_dir=self.assets_path,
            fastq_path='a_fastq_path'
        )
        with open(self.tmp_script) as f:
            assert f.read() == self.writer.script
