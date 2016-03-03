__author__ = 'tcezard'
import logging
from .notification_center import Notification
from analysis_driver.app_logging import logging_default as log_cfg


class LogNotification(Notification):
    """Logs via log_cfg, and also to a notification file with the format'[date time][dataset name] msg'"""
    def __init__(self, dataset, log_file):
        super().__init__(dataset)
        self.ntf_logger = logging.getLogger(self.__class__.__name__)
        handler = logging.FileHandler(filename=log_file, mode='a')
        formatter = logging.Formatter(
            fmt='[%(asctime)s][' + self.dataset.name + '] %(message)s',
            datefmt='%Y-%b-%d %H:%M:%S'
        )
        handler.setFormatter(formatter)
        handler.setLevel(log_cfg.default_level)
        self.ntf_logger.addHandler(handler)

    def start_pipeline(self):
        self._log('info', 'Started pipeline')

    def start_stage(self, stage_name):
        self.dataset.add_stage(stage_name)  # TODO: the dataset should control the notifier
        self._log('info', 'Started stage ' + stage_name)

    def end_stage(self, stage_name, exit_status=0):
        self.dataset.end_stage(stage_name, exit_status)
        if exit_status == 0:
            self._log('info', 'Finished stage ' + stage_name)
        else:
            self._log('error', 'Failed stage ' + stage_name + ' with exit status ' + str(exit_status))

    def end_pipeline(self, exit_status, stacktrace=None):
        self._log('info', 'Finished pipeline with exit status ' + str(exit_status))
        if stacktrace:
            self.ntf_logger.error(self._format_error_message(stacktrace=stacktrace))

    def _log(self, level, msg):
        for l in (self.__logger, self.ntf_logger):
            log_method = l.__getattribute__(level)
            log_method(msg)
