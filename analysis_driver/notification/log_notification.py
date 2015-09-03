__author__ = 'tcezard'
import logging
from analysis_driver.config import logging_default as log_cfg
from analysis_driver.app_logging import AppLogger
from analysis_driver.exceptions import AnalysisDriverError


class LogNotification(AppLogger):
    def __init__(self, log_file):
        self.handler = logging.FileHandler(filename=log_file)
        self.handler.setFormatter(log_cfg.formatter)
        self.handler.setLevel(log_cfg.log_level)

    def start_step(self, step_name, **kwargs):
        self.info('Started {}'.format(step_name))

    def finish_step(self, step_name, run_id, exit_status=0, stop_on_error=False):
        if exit_status == 0:
            self._succeed(step_name)
        else:
            self._fail(step_name, stop_on_error)

    def _succeed(self, step_name, **kwargs):
        self.info('Finished {}'.format(step_name))

    def _fail(self, step_name, stop_on_error):
        self.error('Failed {}'.format(step_name))
        if stop_on_error:
            raise AnalysisDriverError(step_name + ' failed')

    def _check_logger(self):
        if self._logger is None:
            super()._check_logger()
            self._logger.addHandler(self.handler)
