__author__ = 'tcezard'
import pytest
from smtplib import SMTPException
from analysis_driver.notification.notification_center import NotificationCenter
from analysis_driver.notification.email_notification import EmailNotification
from analysis_driver.config import default as cfg
from analysis_driver.exceptions import AnalysisDriverError
from tests.test_analysisdriver import TestAnalysisDriver


class TestNotificationCenter(TestAnalysisDriver):
    def setUp(self):
        self.notification_center = NotificationCenter(cfg, 'test_run')

        email_config = {
            'mailhost': 'renko.ucs.ed.ac.uk',
            'port': 25,
            'reporter_email': 'murray.wham@ed.ac.uk',
            'recipient_emails': ['murray.wham@ed.ac.uk']
        }
        self.email_notification = TestMailRetries('test_run', email_config)

    def test_notification_center(self):
        self.notification_center.start_stage('test stage')

        with pytest.raises(AnalysisDriverError) as e:
            self.notification_center.end_stage('test stage', exit_status=1, stop_on_error=True)
            assert 'test stage failed' in str(e)

    def test_retries(self):
        assert self.email_notification._try_send('this is a test') is True
        assert self.email_notification._try_send('dodgy') is False

        with pytest.raises(AnalysisDriverError) as e:
            self.email_notification._send_mail('dodgy')
            assert 'Failed to send message: dodgy' in str(e)



class TestMailRetries(EmailNotification):

    def _connect_and_send(self, msg):
        """
        Create a test situation where the server connection has failed
        :param email.mime.text.MIMEText msg:
        """
        print('[TestNotificationCenter] Don\'t worry, I\'m not sending any emails.')
        if str(msg).endswith('dodgy'):
            raise SMTPException('Oh noes!')
        else:
            pass
