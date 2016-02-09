__author__ = 'tcezard'
from analysis_driver.app_logging import AppLogger


class _Method:
    # bind a notification center method to its subscribers.
    def __init__(self, send, name):
        self.__send = send
        self.__name = name

    def __call__(self, *args, **kwargs):
        return self.__send(self.__name, *args, **kwargs)


class NotificationCenter(AppLogger):
    """
    Object dispatching notification methods to subscribers
    """
    def __init__(self):
        self.subscribers = []

    def add_subscribers(self, *args):
        """
        e.g: ntf.add_subscribers('a_run_id', (LogNotification, cfg.query('notification', 'log_notification'))
        :param tuple[callable, str, dict] args: A tuple mapping a subscriber class, a run id and a dict
        configuration from analysis_driver.config
        """
        for notifier, run_id, config in args:
            if config:
                self.subscribers.append(notifier(run_id, config))

    def _pass_to_subscriber(self, method_name, *args, **kwargs):
        """
        Take method_name and try to invoke it on each subscriber with *args, **kwargs.
        """
        for subscriber in self.subscribers:
            f = getattr(subscriber, method_name)
            if f and callable(f):
                f(*args, **kwargs)
            else:
                self.debug(
                    'Tried to call nonexistent method %s in %s' % (method_name, subscriber.__class__.__name__)
                )

    def __getattr__(self, name):
        # dispatch a method via getattr, i.e. NotificationCenter.methodname
        return _Method(self._pass_to_subscriber, name)


class Notification(AppLogger):
    def __init__(self, dataset):
        self.dataset = dataset

    def start_pipeline(self):
        pass

    def start_stage(self, stage_name):
        pass

    def end_stage(self, stage_name, exit_status=0):
        pass

    def end_pipeline(self, exit_status, stacktrace=None):
        pass

    @staticmethod
    def _format_error_message(**kwargs):
        msg = ' '.join(('Run failed.', kwargs.get('message', '')))
        if kwargs.get('stacktrace'):
            msg += '\nStack trace below:\n\n'
            msg += kwargs['stacktrace']
        return msg


default = NotificationCenter()
