from unittest.mock import patch
from analysis_driver.dataset import SampleDataset
from tests import TestAnalysisDriver


sample_read_data = 'tests.test_quality_control.qc_tester.SampleDataset._read_data'


class NoCommunicationSample(SampleDataset):
    def _read_data(self):
        pass


class QCTester(TestAnalysisDriver):
    def setUp(self):
        with patch(sample_read_data):
            self.dataset = NoCommunicationSample(
                'test_sample',
                most_recent_proc={
                    'dataset_type': 'sample',
                    'dataset_name': 'test_sample',
                    'proc_id': 'sample_test_sample_now'
                }
            )
