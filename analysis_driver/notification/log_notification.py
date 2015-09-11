__author__ = 'tcezard'
import logging
from analysis_driver.config import logging_default as log_cfg
from analysis_driver.app_logging import AppLogger
from analysis_driver.exceptions import AnalysisDriverError


class LogNotification(AppLogger):
    def __init__(self, run_id, log_file):
        self.run_id = run_id
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
        self.info('Started stage: ' + stage_name)

    def end_stage(self, stage_name, exit_status=0, stop_on_error=False):
        if exit_status == 0:
            self._succeed_stage(stage_name)
        else:
            self._fail_stage(stage_name, exit_status, stop_on_error)

    def end_pipeline(self):
        self.info('Finished pipeline')

    def _succeed_stage(self, stage_name):
        self.info('Finished stage ' + stage_name)

    def _fail_stage(self, stage_name, exit_status, stop_on_error):
        msg = 'Stage ' + stage_name + ' failed with exit status ' + str(exit_status)
        self.error(msg)
        if stop_on_error:
            raise AnalysisDriverError(msg)

    def _check_logger(self):
        if self._logger is None:
            super()._check_logger()  # bind self.logger to the shared handlers...
            self._logger.addHandler(self.handler)  # ... and its own specially-formatted handler
