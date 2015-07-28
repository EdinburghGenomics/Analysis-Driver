__author__ = 'mwham'
from tests.test_writer.test_pbs_writer.test_pbs_writer import TestPBSWriter
from analysis_driver.writer.pbs_writer.fastqc_writer import FastqcWriter


class TestFastqcWriter(TestPBSWriter):

    def setUp(self):
        super().setUp()
        self.writer = FastqcWriter(
            self.tmp_script,
            walltime='24',
            cpus='2',
            mem='1',
            job_name='test',
            log_file='test.log',
            queue='test'
        )

    def test_fastqc(self):
        self.writer._fastqc(self.assets_path)
        for line in ['FASTQ_FILES=`find ' + self.assets_path + ' -name \'*.fastq.gz\'`\n',
                     'fastqc --nogroup -t 8 -q $FASTQ_FILES\n']:
            assert line in self.writer.script

    def test_write(self):
        self.writer.write(input_dir=self.assets_path, run_dir=self.assets_path)
        with open(self.tmp_script) as f:
            assert f.read() == self.writer.script
