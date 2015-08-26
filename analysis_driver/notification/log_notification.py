from analysis_driver.app_logging import AppLogger

__author__ = 'tcezard'



class LogNotification(AppLogger):

    def __init__(self, config):
        pass

    def start_step(self, step_name, **kwargs):
        self.info('Starts {}'.format(step_name))

    def finish_step(self, step_name, **kwargs):
        self.info('Finished {}'.format(step_name))

    def fail_step(self, step_name, **kwargs):
        self.error('Failed {}'.format(step_name))
