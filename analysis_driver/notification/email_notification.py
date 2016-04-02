import os
import smtplib
import jinja2
from email.mime.text import MIMEText
from time import sleep
from .notification_center import Notification
from analysis_driver.exceptions import AnalysisDriverError


class EmailNotification(Notification):
    def __init__(self, dataset, config):
        super().__init__(dataset)
        self.reporter = config['reporter_email']
        self.recipients = config['recipient_emails']
        self.mailhost = config['mailhost']
        self.strict = config.get('strict')
        self.port = config['port']

    def start_pipeline(self):
        self._send_mail('Pipeline started for %s %s ' % (self.dataset.type, self.dataset.name))

    def end_stage(self, stage_name, exit_status=0):
        if exit_status != 0:
            self._send_mail('Stage \'%s\' failed with exit status %s' % (stage_name, exit_status))

    def end_pipeline(self, exit_status=0):
        self._send_mail('Pipeline finished with exit status ' + str(exit_status), bool(exit_status))

    def crash_report(self, exit_status, crash_report):
        self._send_mail(
            'Crash report for %s %s:\n\n%s' % (self.dataset.type, self.dataset.name, crash_report),
            diagnostics=True
        )

    def _send_mail(self, body, diagnostics=False):
        msg = self._prepare_message(body, diagnostics=diagnostics)
        mail_success = self._try_send(msg)
        if not mail_success:
            if self.strict is True:
                raise AnalysisDriverError('Failed to send message: ' + body)
            else:
                self.critical('Failed to send message: ' + body)

    def _try_send(self, msg, retries=3):
        """
        Prepare a MIMEText message from body and diagnostics, and try to send a set number of times.
        :param int retries: Which retry we're currently on
        :return: True if a message is sucessfully sent, otherwise False
        """
        try:
            self._connect_and_send(msg)
            return True
        except (smtplib.SMTPException, TimeoutError) as e:
            retries -= 1
            self.warning('Encountered a %s exception. %s retries remaining', str(e), retries)
            if retries:
                sleep(2)
                return self._try_send(msg, retries)
            else:
                return False

    def _prepare_message(self, body, diagnostics=False):
        """
        Use Jinja to build a MIMEText html-formatted email.
        :param str body: The main body of the email to send
        :param diagnostics: Whether to send diagnostic information (environment variables, configurations)
        """
        content = jinja2.Template(
            open(
                os.path.join(
                    os.path.dirname(os.path.abspath(__file__)),
                    '..', '..', 'etc', 'email_notification.html'
                )
            ).read()
        )
        render_params = {
            'dataset_type': self.dataset.type,
            'run': self.dataset.name,
            'body': self._prepare_string(body, {' ': '&nbsp', '\n': '<br/>'})
        }
        if diagnostics:
            render_params['env_vars'] = self._get_envs('ANALYSISDRIVERCONFIG', 'ANALYSISDRIVERENV')
            # render_params['config'] = self._prepare_string(cfg.report(), {' ': '&nbsp', '\n': '<br/>'})

        msg = MIMEText(content.render(**render_params), 'html')
        msg['Subject'] = 'Analysis Driver %s %s' % (self.dataset.type, self.dataset.name)
        msg['From'] = self.reporter
        msg['To'] = ','.join(self.recipients)

        return msg

    @staticmethod
    def _prepare_string(in_string, charmap):
        for k in charmap:
            in_string = in_string.replace(k, charmap[k])
        return in_string

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
