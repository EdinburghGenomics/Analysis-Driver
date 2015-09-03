__author__ = 'tcezard'
# email modules
import smtplib
import smtpd
from email.mime.text import MIMEText
import subprocess
from analysis_driver.app_logging import AppLogger
from analysis_driver.exceptions import AnalysisDriverError


class EmailNotification(AppLogger):
    def __init__(self, config):
        self.mailhost = config['mailhost']
        # self.reporter_email = config['reporter_email']
        self.recipient_emails = config['recipient_emails']

        # self.server = smtpd.SMTPServer(('localhost', 5000), None)
        # self.connection = smtplib.SMTP(self.mailhost, 465)

    def _send_mail_i(self, subject, body):
        msg = MIMEText(body, 'plain')
        msg['Subject'] = subject
        # msg['From'] = self.reporter_email
        msg['To'] = ','.join(self.recipient_emails)

        p = subprocess.Popen(['mail', '-vt'], stdin=subprocess.PIPE)
        p.communicate(msg.as_bytes())
        exit_status = p.poll()
        if exit_status:
            self.error('Email notification for %s failed with exit status %s' % (subject, exit_status))
        # self.connection.send_message(msg, self.reporter_email, self.recipient_emails)

    def start_step(self, step_name):
        pass

    def finish_step(self, step_name, exit_status=0, stop_on_error=False):
        if exit_status == 0:
            self._succeed(step_name)
        else:
            self._fail(step_name, stop_on_error)

    def _fail(self, step_name, stop_on_error):
        subject = step_name + ' failed'
        body = step_name + ' failed'
        self._send_mail_i(subject, body)
        if stop_on_error:
            raise AnalysisDriverError(subject)

    def _succeed(self, step_name):
        pass
