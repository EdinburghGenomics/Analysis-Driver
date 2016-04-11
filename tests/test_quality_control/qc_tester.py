from unittest.mock import patch
from analysis_driver.dataset_scanner import SampleDataset, MostRecentProc
from tests import TestAnalysisDriver


sample_read_data = 'tests.test_quality_control.qc_tester.SampleDataset._read_data'


class QCTester(TestAnalysisDriver):
    def setUp(self):
        with patch(sample_read_data):
            self.dataset = SampleDataset(
                'test_sample',
                MostRecentProc(
                    'sample',
                    'test_sample',
                    initial_content={
                        'dataset_type': 'sample',
                        'dataset_name': 'test_sample',
                        'proc_id': 'sample_test_sample_now'
                    }
                )
            )
