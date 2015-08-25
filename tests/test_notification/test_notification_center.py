from analysis_driver.notification.notification_center import NotificationCenter, default_notification_center

__author__ = 'tcezard'
from tests.test_analysisdriver import TestAnalysisDriver
from unittest import TestCase

class TestNotificationCenter(TestCase):

    def setUp(self):
        self.notification_center = default_notification_center
    def test_notification_center(self):
        self.notification_center.start_step('test step')