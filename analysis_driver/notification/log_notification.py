__author__ = 'tcezard'
import logging
from .notification_center import Notification
from analysis_driver.config import logging_default as log_cfg


class LogNotification(Notification):
    def __init__(self, dataset, log_file):
        super().__init__(dataset)
        self.handler = logging.FileHandler(filename=log_file, mode='a')
        self.formatter = logging.Formatter(
            fmt='[%(asctime)s][' + self.dataset.name + '] %(message)s',
            datefmt='%Y-%b-%d %H:%M:%S'
        )
        self.handler.setFormatter(self.formatter)
        self.handler.setLevel(log_cfg.log_level)
        # this class will log to the usual places in the usual format, as well as a notification log file in
        # the format '[<date> <time>][dataset name] msg'

    def start_pipeline(self):
        self.info('Started pipeline')

    def start_stage(self, stage_name):
        self.dataset.add_stage(stage_name)
        self.info('Started stage ' + stage_name)

    def end_stage(self, stage_name, exit_status=0):
        self.dataset.remove_stage(stage_name)
        if exit_status == 0:
            self.info('Finished stage ' + stage_name)
        else:
            self.error('Failed stage ' + stage_name + ' with exit status ' + str(exit_status))

    def end_pipeline(self, exit_status, stacktrace=None):
        self.info('Finished pipeline with exit status ' + str(exit_status))
        if stacktrace:
            self.error(self._format_error_message(stacktrace=stacktrace))

    def _check_logger(self):
        """
        Set self._logger as in the superclass, but also bind it to self.handler.
        """
        if self._logger is None:
            super()._check_logger()  # bind self.logger to the shared handlers...
            self._logger.addHandler(self.handler)  # ... and to the differently-formatted self.handler
