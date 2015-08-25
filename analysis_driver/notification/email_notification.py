__author__ = 'tcezard'
import smtplib

# Import the email modules we'll need
from email.mime.text import MIMEText

def _send_mail(mailhost, reporter_email, recipient_emails, subject, body):
    """
    Send out an email with subject and sender if specified.
    :param mailhost: the smtp host as a string
    :param reporter_email: the sender email as a string
    :param notification_emails: a list of emails the message will be sent to
    :param subject: the subject of the message
    :param body: the body of the message
    """
    import smtplib
    msg = MIMEText(body)
    msg['Subject'] = subject
    msg['From'] = reporter_email
    msg['To'] = recipient_emails
    smtpserver = smtplib.SMTP(mailhost, 25)
    smtpserver.send_message(msg, reporter_email, recipient_emails)
    smtpserver.quit()

class EmailNotification:

    def __init__(self, config):
        self._init_with_config(config)

    def _init_with_config(self,config):
        self.mailhost=config['mailhost']
        self.reporter_email=config['reporter_email']
        self.recipient_emails=config['recipient_emails']

    def _send_mail_i(self, subject, body):
        '''convinience method to send email from the Email_notification ofbject'''
        _send_mail(self.mailhost,
                   self.reporter_email,
                   self.recipient_emails,
                   subject, body)

    def start_step(self, step_name, **kwargs):
        subject = "{} has started"
        body = "{} has started"
        self._send_mail_i(subject, body)

    def finish_step(self, step_name, **kwargs):
        subject = "{} has finished"
        body = "{} has finished"
        self._send_mail_i(subject, body)

    def fail_step(self, step_name, **kwargs):
        subject = "{} has failed"
        body = "{} has failed"
        self._send_mail_i(subject, body)

