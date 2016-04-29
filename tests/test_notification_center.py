import pytest
from unittest.mock import Mock, patch
from smtplib import SMTPException
from analysis_driver.dataset_scanner import NoCommunicationDataset
from analysis_driver.notification.notification_center import NotificationCenter
from analysis_driver.notification import EmailNotification, LogNotification, AsanaNotification
from analysis_driver.config import default as cfg
from analysis_driver.exceptions import AnalysisDriverError
from tests.test_analysisdriver import TestAnalysisDriver


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
        dataset = NoCommunicationDataset('test_run_id')
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
        self.email_ntf = EmailNotification(dataset, email_config)

    @patch('smtplib.SMTP', new=FakeSMTP)
    @patch('analysis_driver.notification.email_notification.sleep')
    def test_retries(self, mocked_time):
        assert self.email_ntf._try_send(self.email_ntf._prepare_message('this is a test')) is True
        assert self.email_ntf._try_send(self.email_ntf._prepare_message('dodgy')) is False

        with pytest.raises(AnalysisDriverError) as e:
            self.email_ntf._send_mail('dodgy')
            assert 'Failed to send message: dodgy' in str(e)


class TestAsanaNotification(TestAnalysisDriver):
    def setUp(self):
        self.ntf = AsanaNotification(
            NoCommunicationDataset('test_dataset'),
            {'access_token': 'an_access_token', 'workspace_id': 1337, 'project_id': 1338}
        )
        self.ntf.client = Mock(
            tasks=Mock(
                find_all=Mock(return_value=[{'name': 'this'}]),
                create_in_workspace=Mock(return_value={'id': 1337}),
                find_by_id=Mock(return_value={'name': 'this', 'id': 1337})
            )
        )

    def test_task(self):
        assert self.ntf._task is None
        assert self.ntf.task == {'id': 1337, 'name': 'this'}
        self.ntf.client.tasks.find_by_id.assert_called_with(1337)

    def test_add_comment(self):
        self.ntf._add_comment('a comment')
        self.ntf.client.tasks.add_comment.assert_called_with(1337, text='a comment')

    def test_get_entity(self):
        collection = [{'name': 'this'}, {'name': 'that'}]
        assert self.ntf._get_entity(collection, 'that') == {'name': 'that'}
        assert self.ntf._get_entity(collection, 'other') is None

    def test_create_task(self):
        assert self.ntf._create_task() == {'id': 1337}
        self.ntf.client.tasks.create_in_workspace.assert_called_with(1337, self.ntf.task_template)
