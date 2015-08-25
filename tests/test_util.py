__author__ = 'mwham'
from tests.test_analysisdriver import TestAnalysisDriver
from analysis_driver.util import fastq_handler
from analysis_driver.app_logging import AppLogger
import os.path

helper = TestAnalysisDriver()


def test_setup_bcbio_run():
    print('Currently untestable')
    assert True


class TestLogger(TestAnalysisDriver):
    def setUp(self):
        self.logger = AppLogger()

    def test(self):
        assert True  # Not much to test here


class TestFastqHandler(TestAnalysisDriver):

    def test_find_fastqs(self):
        fastqs = fastq_handler.find_fastqs(self.fastq_path, '10015AT')
        for file_name in ['this.fastq.gz', 'that.fastq.gz', 'other.fastq.gz']:
            assert os.path.join(
                self.fastq_path, '10015AT', '10015ATA0001L05', file_name
            ) in fastqs['10015ATA0001L05']

    def test_flatten_fastqs(self):
        fastqs = fastq_handler.flatten_fastqs(self.fastq_path, ['10015AT'])
        for file_name in ['this.fastq.gz', 'that.fastq.gz', 'other.fastq.gz']:
            assert os.path.join(
                self.fastq_path, '10015AT', '10015ATA0001L05', file_name
            ) in fastqs
