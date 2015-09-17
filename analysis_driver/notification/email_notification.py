__author__ = 'tcezard'
import smtplib
from email.mime.text import MIMEText
from time import sleep
from analysis_driver.app_logging import AppLogger
from analysis_driver.exceptions import AnalysisDriverError


class EmailNotification(AppLogger):
    def __init__(self, run_id, cfg):
        self.run_id = run_id
        self.reporter = cfg['reporter_email']
        self.recipients = cfg['recipient_emails']
        self.mailhost = cfg['mailhost']
        self.port = cfg['port']

    def start_pipeline(self):
        self._send_mail('Pipeline started for run ' + self.run_id)

    def start_stage(self, stage_name):
        pass

    def end_stage(self, stage_name, exit_status=0, stop_on_error=False):
        if exit_status == 0:
            pass
        else:
            self._fail_stage(stage_name, exit_status, stop_on_error)

    def end_pipeline(self):
        self._send_mail('Pipeline finished for run ' + self.run_id)

    def _fail_stage(self, stage_name, exit_status, stop_on_error):
        msg = ('%s failed for run %s with exit status %s' % (stage_name, self.run_id, exit_status))
        self._send_mail(msg)
        if stop_on_error:
            raise AnalysisDriverError(msg)

    def _send_mail(self, body):
        mail_success = self._try_send(body)
        if not mail_success:
            self.critical('Failed to send message: ' + body, error_class=AnalysisDriverError)

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
        msg = MIMEText(body, 'plain')
        msg['Subject'] = 'Run ' + self.run_id
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
