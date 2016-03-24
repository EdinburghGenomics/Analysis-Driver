from analysis_driver.exceptions import AnalysisDriverError

__author__ = 'tcezard'
from tests.test_analysisdriver import TestAnalysisDriver
from unittest.mock import patch

from bin.deliver_reviewed_data import DataDelivery


patched_deliverable_project = patch(
    'analysis_driver.rest_communication.get_documents',
    return_value=[
        {
            'sample_id': 'deliverable_sample',
            'project_id': 'test_project',
            'user_sample_id': 'user_s_id',
            'useable': 'yes',
            'analysis_driver_procs': [
                {
                    'proc_id': 'most_recent_proc',
                    '_created': '06_03_2016_12:00:00',
                    'status': 'finished'
                },
                {
                    'proc_id': 'old_recent_proc',
                    '_created': '01_02_2016_12:00:00',
                    'status': 'aborted'
                }
            ],
            'bam_file_reads': 1,
            'run_elements' : [
                {

                }
            ]
        }
    ]
)
patched_error_project = patch(
    'analysis_driver.rest_communication.get_documents',
    return_value=[
        {
            'sample_id': 'deliverable_sample',
            'project_id': 'test_project',
            'useable': 'yes',
            'analysis_driver_procs': [
                {
                    'proc_id': 'most_recent_proc',
                    '_created': '06_03_2016_12:00:00',
                    'status': 'aborted'
                },
                {
                    'proc_id': 'old_recent_proc',
                    '_created': '01_02_2016_12:00:00',
                    'status': 'finished'
                }
            ]
        }
    ]
)


class TestDataDelivery(TestAnalysisDriver):

    def setUp(self):
        self.delivery_dry = DataDelivery(dry_run=True)

    def test_get_deliverable_projects_samples(self):
        with patched_deliverable_project as mocked_get_doc:
            project_to_samples = self.delivery_dry.get_deliverable_projects_samples(project_id='test_project')
            assert dict(project_to_samples) == {'test_project': ['deliverable_sample']}

        with patched_error_project as mocked_get_doc:
            self.assertRaises(AnalysisDriverError, self.delivery_dry.get_deliverable_projects_samples)



    def test_summarise_metrics_per_sample(self):
        with patched_deliverable_project as mocked_get_doc:
            self.delivery_dry.get_deliverable_projects_samples(project_id='test_project')
            expected_header = ['Project', 'Sample Id', 'User sample id', 'Read pair sequenced',
                               'Yield', 'Yield Q30', 'Nb reads in bam', 'mapping rate', 'properly mapped reads rate',
                               'duplicate rate', 'Mean coverage', 'Callable bases rate', 'Delivery folder']
            expected_lines=[
                'test_project\tdeliverable_sample\tuser_s_id\t0\t0.0\t0.0\t1\t0.0\t0.0\t0.0\t0.0\t0\tdate_delivery'
            ]
            with patch('analysis_driver.clarity.get_species_from_sample', return_value="Homo sapiens"):
                header, lines = self.delivery_dry.summarise_metrics_per_sample(project_id='test_project', delivery_folder='date_delivery')
                assert header == expected_header
                assert lines == expected_lines
