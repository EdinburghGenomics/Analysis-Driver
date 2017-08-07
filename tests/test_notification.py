import pytest
from unittest.mock import patch, Mock
from analysis_driver import notification
from tests.test_analysisdriver import TestAnalysisDriver
from pyclarity_lims.lims import Lims
from pyclarity_lims.entities import Step

class TestNotificationCentre(TestAnalysisDriver):

    def setUp(self):
        self.n = notification.NotificationCentre('name')

    @patch('egcg_core.notifications.NotificationCentre.notify')
    def test_start_pipeline(self, mocked_notify):
        self.n.start_pipeline()
        mocked_notify.assert_called_with('Started pipeline', ('log', 'email'))

    @patch('egcg_core.notifications.NotificationCentre.notify')
    def test_start_stage(self, mocked_notify):
        self.n.start_stage('test_stage')
        mocked_notify.assert_called_with('Started stage test_stage', ('log',))

    @patch('egcg_core.notifications.NotificationCentre.notify')
    def test_end_stage(self, mocked_notify):
        self.n.end_stage('test_stage')
        mocked_notify.assert_called_with("Stage 'test_stage' finished with exit status 0", ['log'])

    @patch('egcg_core.notifications.NotificationCentre.notify')
    def test_end_pipeline(self, mocked_notify):
        self.n.end_pipeline(0)
        mocked_notify.assert_called_with('Pipeline finished with exit status 0', ['log', 'email'])

    @patch('egcg_core.notifications.NotificationCentre.notify_all')
    def test_crash_report(self, mocked_notify_all):
        self.n.crash_report('stacktrace')
        mocked_notify_all.assert_called_with('stacktrace')


class TestLimsNotification(TestAnalysisDriver):
    def setUp(self):
        with patch('egcg_core.clarity.get_sample', return_value=Mock(artifact='LimsArtifact')),\
             patch('egcg_core.clarity.get_workflow_stage', return_value=Mock(uri='stage_uri', step='stage_step')):
            self.n = notification.LimsNotification('test_name')

    @patch.object(Lims, 'route_artifacts')
    def test_route_sample_to_data_processing(self, mocked_route):
        self.n.route_sample_to_data_processing()
        mocked_route.assert_called_with(['LimsArtifact'], stage_uri='stage_uri')

    @patch('egcg_core.clarity.connection', return_value='LIMSconnection')
    @patch.object(Step, 'create')
    def test_create_step(self, mocked_create_step, mocked_lims):
        self.n.create_step()
        mocked_create_step.assert_called_with('LIMSconnection', inputs='LimsArtifact', protocol_step='stage_step')
        mocked_create_step.advance.assert_called_once()

    def test_assign_next_and_advance_step(self):
        self.n.step = Mock(current_state='Assign Next Steps', actions=Mock(next_actions=[{}]))
        self.n.assign_next_and_advance_step()
        self.n.step.actions.put.assert_called_once()
        self.n.step.advance.assert_called_once()
        assert self.n.step.actions.next_actions[0]['action'] == 'nextstep'
        assert self.n.step.actions.next_actions[0]['step-uri'] == 'stage_uri'
        self.n.step = Mock(current_state='Incorrect State', actions=Mock(next_actions=[{}]))
        with pytest.raises(AssertionError):
            self.n.assign_next_and_advance_step()

    def test_remove_sample_from_workflow(self):
        self.n.step = Mock(current_state='Assign Next Steps', actions=Mock(next_actions=[{}]))
        self.n.remove_sample_from_workflow()
        self.n.step.actions.put.assert_called_once()
        self.n.step.advance.assert_called_once()
        assert self.n.step.actions.next_actions[0]['action'] == 'remove'
        self.n.step = Mock(current_state='Incorrect State', actions=Mock(next_actions=[{}]))
        with pytest.raises(AssertionError):
            self.n.remove_sample_from_workflow()

    @patch('analysis_driver.notification.LimsNotification.route_sample_to_data_processing')
    @patch('analysis_driver.notification.LimsNotification.create_step')
    def test_start_sample_pipeline(self, mocked_create_step, mocked_route):
        self.n.start_sample_pipeline()
        mocked_route.assert_called_once()
        mocked_create_step.assert_called_once()
