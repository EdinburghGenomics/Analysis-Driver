import os
from analysis_driver.dataset_scanner import RunDataset

__author__ = 'tcezard'
import pytest
import sys
from smtplib import SMTPException
from analysis_driver.notification.notification_center import NotificationCenter
from analysis_driver.notification import EmailNotification, LogNotification
from analysis_driver.config import default as cfg
from analysis_driver.exceptions import AnalysisDriverError
from tests.test_analysisdriver import TestAnalysisDriver

if not sys.argv:
    print('Usage: python test_notification_center.py <mailhost> <port> <reporter_email> <recipient_emails>')
    sys.exit(0)


class TestNotificationCenter(TestAnalysisDriver):
    def setUp(self):
        base_dir = os.path.join(self.assets_path, 'dataset_scanner')
        dataset = RunDataset(
            name='test_run_id',
            path=os.path.join(base_dir, 'that'),
            lock_file_dir=base_dir,
            use_int_dir=False
        )
        self.notification_center = NotificationCenter()
        self.notification_center.add_subscribers(
            (LogNotification, dataset, cfg.query('notification', 'log_notification')),
            (TestEmailNotification, dataset, cfg.query('notification', 'email_notification'))
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
        self.email_notification = TestEmailNotification(dataset, email_config)
        if cfg.query('notification', 'email_notification'):
            cfg.content['notification']['email_notification']['strict'] = True

    def test_retries(self):
        assert self.email_notification._try_send('this is a test', diagnostics=False) is True
        assert self.email_notification._try_send('dodgy', diagnostics=False) is False

        with pytest.raises(AnalysisDriverError) as e:
            self.email_notification._send_mail('dodgy')
            assert 'Failed to send message: dodgy' in str(e)


class TestEmailNotification(EmailNotification):

    def _connect_and_send(self, msg):
        """
        Create a test situation where the server connection has failed
        :param email.mime.text.MIMEText msg:
        """
        print('[TestNotificationCenter] Don\'t worry, I\'m not sending any emails.')
        if 'dodgy' in str(msg):
            raise SMTPException('Oh noes!')
        else:
            pass
