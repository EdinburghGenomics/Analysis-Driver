__author__ = 'tcezard'
import os
import pytest
from unittest.mock import Mock, patch
from smtplib import SMTPException
from analysis_driver.dataset_scanner import RunDataset
from analysis_driver.notification.notification_center import NotificationCenter
from analysis_driver.notification import EmailNotification, LogNotification
from analysis_driver.config import default as cfg
from analysis_driver.exceptions import AnalysisDriverError
from tests.test_analysisdriver import TestAnalysisDriver
from tests.test_dataset_scanner import patched_request


class FakeSMTP(Mock):
    def __init__(self, host, port):
        super().__init__()
        self.mailhost, self.port = host, port

    @staticmethod
    def send_message(msg, reporter, recipients):
        if 'dodgy' in str(msg):
            raise SMTPException('Oh noes!')
        else:
            pass

    @staticmethod
    def quit():
        pass


class TestNotificationCenter(TestAnalysisDriver):
    def setUp(self):
        base_dir = os.path.join(self.assets_path, 'dataset_scanner')
        with patched_request:
            dataset = RunDataset(
                name='test_run_id',
                path=os.path.join(base_dir, 'that'),
                use_int_dir=False
            )
        self.notification_center = NotificationCenter()
        self.notification_center.add_subscribers(
            (LogNotification, dataset, cfg.query('notification', 'log_notification')),
            (EmailNotification, dataset, cfg.query('notification', 'email_notification'))
        )

        email_config = {
            'strict': True,
            'mailhost': 'a_mailhost',
            'port': 0,
            'reporter_email': 'a_reporter',
            'recipient_emails': ['a_recipient', 'another_recipient']
        }
        if cfg.query('notification', 'email_notification'):
            email_config = cfg.query('notification', 'email_notification')
        
        print(self.notification_center.subscribers)
        self.email_notification = EmailNotification(dataset, email_config)

    @patch('smtplib.SMTP', new=FakeSMTP)
    def test_retries(self):
        assert self.email_notification._try_send('this is a test', diagnostics=False) is True
        assert self.email_notification._try_send('dodgy', diagnostics=False) is False

        with pytest.raises(AnalysisDriverError) as e:
            self.email_notification._send_mail('dodgy')
            assert 'Failed to send message: dodgy' in str(e)