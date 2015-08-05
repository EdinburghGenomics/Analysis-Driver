__author__ = 'mwham'
from tests.test_writer.test_pbs_writer.test_pbs_writer import TestPBSWriter
from analysis_driver.writer.pbs_writer.fastqc_writer import FastqcWriter
from analysis_driver.util.fastq_handler import flatten_fastqs


class TestFastqcWriter(TestPBSWriter):

    def setUp(self):
        super().setUp()
        self.writer = FastqcWriter(
            self.tmp_script,
            job_name='test',
            fastqs=flatten_fastqs(self.fastq_path, ['10015AT']),
            walltime='24',
            cpus='2',
            mem='1',
            log_file='test.log',
            queue='test'
        )

    def test_fastqc(self):
        self.writer._fastqc(self.writer.fastqs)
        for line in [
            '#PBS -J 1-3\n',
            'case $PBS_ARRAY_INDEX in\n',
            '1) fastqc --nogroup -q /Users/mwham/workspace/EdGen_Analysis_Driver/Applications/Analysis-Driver/tests/assets/fastqs/10015AT/10015ATA0001L05/other.fastq.gz\n',
            ';;\n',
            '2) fastqc --nogroup -q /Users/mwham/workspace/EdGen_Analysis_Driver/Applications/Analysis-Driver/tests/assets/fastqs/10015AT/10015ATA0001L05/that.fastq.gz\n',
            '3) fastqc --nogroup -q /Users/mwham/workspace/EdGen_Analysis_Driver/Applications/Analysis-Driver/tests/assets/fastqs/10015AT/10015ATA0001L05/this.fastq.gz\n',
            '*) echo "Unexpected PBS_ARRAY_INDEX: $PBS_ARRAY_INDEX"\n',
            'esac\n'
        ]:
            assert line in self.writer.script

    def test_write(self):
        self.writer.write()
        with open(self.tmp_script) as f:
            assert f.read() == self.writer.script
