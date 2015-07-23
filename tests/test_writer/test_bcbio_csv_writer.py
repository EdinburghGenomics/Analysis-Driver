__author__ = 'mwham'
from tests import TestAnalysisDriver
from analysis_driver.writer import BCBioCSVWriter
from analysis_driver.reader.sample_sheet import SampleSheet
import os


class TestBCBioCSVWriter(TestAnalysisDriver):
    def setUp(self):
        sample_sheet = SampleSheet(self.assets_path)
        self.tmp_run_dir = os.path.join(os.path.dirname(__file__), 'test_run_dir')
        if not os.path.exists(self.tmp_run_dir):
            os.makedirs(self.tmp_run_dir)
        self.writer = BCBioCSVWriter(self.fastq_path, self.tmp_run_dir, sample_sheet)

    def test_write(self):
        self.writer.write()
        with open(os.path.join(self.tmp_run_dir, 'samples.csv')) as f:
            print(f)
            assert True
            # TODO: assert that the file has correct content

    def test_find_fastqs(self):
        fastqs = self.writer._find_fastqs(self.fastq_path, '10015AT')
        for file_name in ['this.fastq.gz', 'that.fastq.gz', 'other.fastq.gz']:
            assert os.path.join(self.fastq_path, '10015AT', '10015ATA0001L05', file_name) in fastqs
