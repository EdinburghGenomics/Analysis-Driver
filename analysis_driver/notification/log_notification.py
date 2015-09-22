__author__ = 'tcezard'
import logging
from .notification_center import Notification
from analysis_driver.config import logging_default as log_cfg


class LogNotification(Notification):
    def __init__(self, run_id, log_file):
        super().__init__(run_id)
        self.handler = logging.FileHandler(filename=log_file)
        self.formatter = logging.Formatter(
            fmt='[%(asctime)s][' + self.run_id + '] %(message)s',
            datefmt='%Y-%b-%d %H:%M:%S'
        )
        self.handler.setFormatter(self.formatter)
        self.handler.setLevel(log_cfg.log_level)

    def start_pipeline(self):
        self.info('Started pipeline')

    def start_stage(self, stage_name):
        self.info('Started stage ' + stage_name)

    def end_stage(self, stage_name, exit_status=0):
        if exit_status == 0:
            self.info('Finished stage ' + stage_name)
        else:
            self.error('Failed stage ' + stage_name + ' with exit status ' + str(exit_status))

    def end_pipeline(self, exit_status, stacktrace=None):
        self.info('Finished pipeline with exit status ' + str(exit_status))
        if stacktrace:
            self.error(self._format_error_message(stacktrace=stacktrace))

    def _check_logger(self):
        if self._logger is None:
            super()._check_logger()  # bind self.logger to the shared handlers...
            self._logger.addHandler(self.handler)  # ... and its own specially-formatted handler
