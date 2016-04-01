from unittest.mock import patch
from analysis_driver.dataset_scanner import SampleDataset
from tests import TestAnalysisDriver


sample_read_data = 'tests.test_quality_control.qc_tester.SampleDataset._read_data'
sample_most_recent_proc = 'tests.test_quality_control.qc_tester.SampleDataset._most_recent_proc'


class QCTester(TestAnalysisDriver):
    def setUp(self):
        with patch(sample_read_data), patch(sample_most_recent_proc):
            self.dataset = SampleDataset('test_sample')
