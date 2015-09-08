__author__ = 'tcezard'
import smtplib
from email.mime.text import MIMEText
from analysis_driver.app_logging import AppLogger
from analysis_driver.exceptions import AnalysisDriverError


class EmailNotification(AppLogger):
    def __init__(self, cfg):
        self.reporter = cfg['reporter_email']
        self.recipients = cfg['recipient_emails']
        self.mailhost = cfg['mailhost']
        self.port = cfg['port']

    def start_pipeline(self, run_id):
        self._send_mail('Run ' + run_id, 'Pipeline started for run ' + run_id)

    def start_stage(self, stage_name):
        pass

    def end_stage(self, stage_name, run_id, exit_status=0, stop_on_error=False):
        if exit_status == 0:
            pass
        else:
            self._fail_stage(stage_name, run_id, stop_on_error)

    def end_pipeline(self, run_id):
        self._send_mail('Run ' + run_id, 'Pipeline finished for run ' + run_id)

    def _fail_stage(self, stage_name, run_id, stop_on_error):
        self._send_mail('Run ' + run_id, stage_name + ' failed for run ' + run_id)
        if stop_on_error:
            raise AnalysisDriverError(stage_name + ' failed')

    def _send_mail(self, subject, body):
        connection = smtplib.SMTP(self.mailhost, self.port)
        msg = MIMEText(body, 'plain')
        msg['Subject'] = subject
        msg['From'] = self.reporter
        msg['To'] = ','.join(self.recipients)

        connection.send_message(
            msg,
            self.reporter,
            self.recipients
        )
        connection.quit()

