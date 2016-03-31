import logging
from .notification_center import Notification
from analysis_driver.app_logging import logging_default as log_cfg


class LogNotification(Notification):
    """Logs via log_cfg, and also to a notification file with the format'[date time][dataset_name] msg'"""
    def __init__(self, dataset, log_file):
        super().__init__(dataset)
        handler = logging.FileHandler(filename=log_file, mode='a')
        handler.setFormatter(
            logging.Formatter(
                fmt='[%(asctime)s][' + self.dataset.name + '] %(message)s',
                datefmt='%Y-%b-%d %H:%M:%S'
            )
        )
        handler.setLevel(log_cfg.log_level)
        self._logger.addHandler(handler)

    def start_pipeline(self):
        self.info('Started pipeline')

    def start_stage(self, stage_name):
        self.info('Started stage ' + stage_name)

    def end_stage(self, stage_name, exit_status=0):
        self.info('Finished stage ' + stage_name + ' with exit status ' + str(exit_status))

    def end_pipeline(self, exit_status=0):
        self.info('Finished pipeline with exit status ' + str(exit_status))
