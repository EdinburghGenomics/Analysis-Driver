from unittest.mock import Mock
from analysis_driver.dataset import SampleDataset, RunDataset
from tests import TestAnalysisDriver


class QCTester(TestAnalysisDriver):
    def setUp(self):
        self.dataset = SampleDataset(
            'test_sample',
            most_recent_proc={
                'dataset_type': 'sample',
                'dataset_name': 'test_sample',
                'proc_id': 'sample_test_sample_now'
            }
        )

        self.run_dataset = RunDataset(
            'test_run',
            most_recent_proc={
                'dataset_type': 'run',
                'dataset_name': 'test_run',
                'proc_id': 'run_test_run_now'
            }
        )

        self.dataset.most_recent_proc = Mock()
        self.dataset._ntf = Mock()
        self.run_dataset.most_recent_proc = Mock()
        self.run_dataset._ntf = Mock()
