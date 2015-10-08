__author__ = 'mwham'
from tests.test_analysisdriver import TestAnalysisDriver
from analysis_driver import util
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
        fastqs = util.fastq_handler.find_fastqs(self.fastq_path, '10015AT', '10015AT0001')
        for file_name in ['10015AT0001_merged_R1.fastq.gz', '10015AT0001_merged_R2.fastq.gz']:
            assert os.path.join(
                self.fastq_path, '10015AT', '10015AT0001', file_name
            ) in fastqs

    def test_find_all_fastqs(self):
        fastqs = util.fastq_handler.find_all_fastqs(self.fastq_path)
        for file_name in ['10015AT0001_merged_R1.fastq.gz', '10015AT0001_merged_R2.fastq.gz']:
            assert os.path.join(
                self.fastq_path, '10015AT', '10015AT0001', file_name
            ) in fastqs


def test_transfer_output_files():
    sample_id = '10015AT0001'
    destination = os.path.join(helper.data_output, 'output_data')
    if not os.path.isdir(destination):
        os.mkdir(destination)
    for f in os.listdir(destination):
        os.remove(os.path.join(destination, f))
    assert not os.listdir(destination)

    source_path_mapping = {
        'vcf': os.path.join(helper.data_output, 'samples_10015AT0001-merged', 'final'),
        'bam': os.path.join(helper.data_output, 'samples_10015AT0001-merged', 'final'),
        'fastq': os.path.join(helper.data_output, 'merged_fastqs')
    }

    util.transfer_output_files(sample_id, destination, source_path_mapping)

    output_files = os.listdir(destination)
    output_files.sort()

    expected_outputs = [
        '10015AT0001.bam',
        '10015AT0001.bam.bai',
        '10015AT0001.bam.bai.md5',
        '10015AT0001.bam.md5',
        '10015AT0001.g.vcf.gz',
        '10015AT0001.g.vcf.gz.md5',
        '10015AT0001.g.vcf.gz.tbi',
        '10015AT0001.g.vcf.gz.tbi.md5',
        '10015AT0001_R1.fastq.gz',
        '10015AT0001_R1.fastq.gz.md5',
        '10015AT0001_R2.fastq.gz',
        '10015AT0001_R2.fastq.gz.md5'
        ]
    print(output_files)
    print(expected_outputs)
    assert output_files == expected_outputs
