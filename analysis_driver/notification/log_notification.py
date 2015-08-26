__author__ = 'tcezard'
import logging
from analysis_driver.config import logging_default as log_cfg
from analysis_driver.app_logging import AppLogger


class LogNotification(AppLogger):
    def __init__(self, log_file):
        self.handler = logging.FileHandler(filename=log_file)
        self.handler.setFormatter(log_cfg.formatter)
        self.handler.setLevel(log_cfg.log_level)

    def start_step(self, step_name, **kwargs):
        self.info('Started {}'.format(step_name))

    def finish_step(self, step_name, **kwargs):
        self.info('Finished {}'.format(step_name))

    def fail_step(self, step_name, **kwargs):
        self.error('Failed {}'.format(step_name))

    def _check_logger(self):
        if self._logger is None:
            super()._check_logger()
            self._logger.addHandler(self.handler)
