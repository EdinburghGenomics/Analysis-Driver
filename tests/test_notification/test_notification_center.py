__author__ = 'tcezard'
import pytest
from analysis_driver.notification.notification_center import default_notification_center
from analysis_driver.exceptions import AnalysisDriverError
from tests.test_analysisdriver import TestAnalysisDriver


class TestNotificationCenter(TestAnalysisDriver):

    def setUp(self):
        self.notification_center = default_notification_center

    def test_notification_center(self):
        self.notification_center.start_stage('test stage')

        with pytest.raises(AnalysisDriverError) as e:
            self.notification_center.end_stage('test stage', 'a_run_id', exit_status=1, stop_on_error=True)
            assert 'test stage failed' in str(e)
