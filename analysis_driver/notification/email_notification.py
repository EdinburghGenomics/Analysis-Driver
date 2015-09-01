__author__ = 'tcezard'
import smtplib
# Import the email modules we'll need
from email.mime.text import MIMEText
from analysis_driver.exceptions import AnalysisDriverError


def _send_mail(mailhost, reporter_email, recipient_emails, subject, body):
    """
    Send out an email with subject and sender if specified.
    :param mailhost: the smtp host as a string
    :param reporter_email: the sender email as a string
    :param recipient_emails: a list of emails the message will be sent to
    :param subject: the subject of the message
    :param body: the body of the message
    """
    msg = MIMEText(body)
    msg['Subject'] = subject
    msg['From'] = reporter_email
    msg['To'] = recipient_emails
    smtpserver = smtplib.SMTP(mailhost, 25)
    smtpserver.send_message(msg, reporter_email, recipient_emails)
    smtpserver.quit()


class EmailNotification:
    def __init__(self, config):
        self.mailhost = config['mailhost']
        self.reporter_email = config['reporter_email']
        self.recipient_emails = config['recipient_emails']

    def _send_mail_i(self, subject, body):
        """
        Convenience method allowing EmailNotification to call _send_mail
        """
        _send_mail(
            self.mailhost,
            self.reporter_email,
            self.recipient_emails,
            subject,
            body
        )

    def start_step(self, step_name):
        subject = step_name + ' started'
        body = step_name + ' started'
        self._send_mail_i(subject, body)

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
        subject = step_name + ' finished'
        body = step_name + ' finished'
        self._send_mail_i(subject, body)
