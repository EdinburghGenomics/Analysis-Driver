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

    def start_pipeline(self, *args):
        self.info('Started pipeline')

    def start_stage(self, stage_name):
        self.info('Started stage: ' + stage_name)

    def end_stage(self, stage_name, run_id, exit_status=0, stop_on_error=False):
        if exit_status == 0:
            self._succeed_stage(stage_name)
        else:
            self._fail_stage(stage_name, stop_on_error)

    def end_pipeline(self, *args):
        self.info('Finished pipeline')

    def _succeed_stage(self, stage_name):
        self.info('Finished stage ' + stage_name)

    def _fail_stage(self, stage_name, stop_on_error):
        self.error('Failed stage ' + stage_name)
        if stop_on_error:
            raise AnalysisDriverError(stage_name + ' failed')

    def _check_logger(self):
        if self._logger is None:
            super()._check_logger()
            self._logger.addHandler(self.handler)
