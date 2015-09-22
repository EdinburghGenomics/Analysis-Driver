__author__ = 'tcezard'
import pytest
from smtplib import SMTPException
from analysis_driver.notification.notification_center import NotificationCenter
from analysis_driver.notification import EmailNotification, LogNotification
from analysis_driver.config import default as cfg
from analysis_driver.exceptions import AnalysisDriverError
from tests.test_analysisdriver import TestAnalysisDriver


class TestNotificationCenter(TestAnalysisDriver):
    def setUp(self):
        self.notification_center = NotificationCenter()
        self.notification_center.add_subscribers(
            (LogNotification, 'test_run_id', cfg.query('notification', 'log_notification')),
            (TestEmailNotification, 'test_run_id', cfg.query('notification', 'email_notification'))
        )

        email_config = {
            'mailhost': 'renko.ucs.ed.ac.uk',
            'port': 25,
            'reporter_email': 'murray.wham@ed.ac.uk',
            'recipient_emails': ['murray.wham@ed.ac.uk']
        }
        print(self.notification_center.subscribers)
        self.email_notification = TestEmailNotification('test_run', email_config)

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
