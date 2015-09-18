__author__ = 'tcezard'
import smtplib
from email.mime.text import MIMEText
import jinja2
import os.path
from time import sleep
from analysis_driver.config import default as cfg
from .notification_center import Notification
from analysis_driver.exceptions import AnalysisDriverError


class EmailNotification(Notification):
    def __init__(self, run_id, config):
        super().__init__(run_id)
        self.reporter = config['reporter_email']
        self.recipients = config['recipient_emails']
        self.mailhost = config['mailhost']
        self.port = config['port']

    def start_pipeline(self):
        self._send_mail('Pipeline started for run ' + self.run_id)

    def end_stage(self, stage_name, exit_status=0):
        if exit_status == 0:
            pass
        else:
            self._send_mail(
                'Stage \'%s\' failed with exit status %s' % (stage_name, exit_status)
            )

    def end_pipeline(self):
        self._send_mail('Pipeline finished for run ' + self.run_id)

    def fail_pipeline(self, message='', **kwargs):
        self._send_mail(self._format_error_message(message, kwargs.get('stacktrace')))

    def _send_mail(self, body):
        mail_success = self._try_send(body)
        if not mail_success:
            raise AnalysisDriverError('Failed to send message: ' + body)

    def _try_send(self, body, retries=1):
        msg = self._prepare_message(body)
        try:
            self._connect_and_send(msg)
            return True
        except smtplib.SMTPException as e:
            self.warn('Encountered a ' + str(e) + ' exception. Retry number ' + str(retries))
            retries += 1
            if retries <= 3:
                sleep(2)
                return self._try_send(body, retries)
            else:
                return False

    def _prepare_message(self, body):
        content = jinja2.Template(
            open(
                os.path.join(
                    os.path.dirname(os.path.abspath(__file__)),
                    '..', '..', 'etc', 'email_notification.html'
                )
            ).read()
        )
        msg = MIMEText(
            content.render(
                run=self.run_id,
                body=body,
                env_vars=self._get_envs('ANALYSISDRIVERCONFIG', 'ANALYSISDRIVERENV'),
                run_config=cfg.report(space='&nbsp')
            ),
            'html'
        )
        msg['Subject'] = 'Analysis Driver run ' + self.run_id
        msg['From'] = self.reporter
        msg['To'] = ','.join(self.recipients)

        return msg

    def _connect_and_send(self, msg):
        connection = smtplib.SMTP(self.mailhost, self.port)
        connection.send_message(
            msg,
            self.reporter,
            self.recipients
        )
        connection.quit()

    @staticmethod
    def _get_envs(*envs):
        return ((e, os.getenv(e)) for e in envs)
