__author__ = 'tcezard'
from analysis_driver.app_logging import AppLogger


class _Method:
    # some magic to bind an notification center method to its subscribers.
    def __init__(self, send, name):
        self.__send = send
        self.__name = name

    def __call__(self, *args, **kwargs):
        return self.__send(self.__name, *args, **kwargs)


class NotificationCenter(AppLogger):
    """
    Object dispatching notification to subscribers.
    The suscribers are defined in a config file
    """
    def __init__(self):
        self.subscribers = []

    def add_subscribers(self, run_id, *subscribers):
        """
        e.g: ntf.add_subscribers('a_run_id', (LogNotification, cfg.query('notification', 'log_notification'))
        :param str run_id: A run id to assign to the subscriber
        :param tuple[callable, dict] subscribers: This should be a tuple containing the ClassName and
        configuration
        """
        for notifier, config in subscribers:
            if config:
                self.subscribers.append(notifier(run_id, config))

    def _pass_to_subscriber(self, function_name, *args, **kwargs):
        for subscriber in self.subscribers:
            f = getattr(subscriber, function_name)
            if f and callable(f):
                f(*args, **kwargs)
            else:
                self.debug(
                    'Tried to call nonexistent method %s in class %s' % (
                        function_name, subscriber.__class__.__name__
                    )
                )

    def __getattr__(self, name):
        # magic method dispatcher
        return _Method(self._pass_to_subscriber, name)


class Notification(AppLogger):
    def __init__(self, run_id):
        self.run_id = run_id

    def start_pipeline(self):
        pass

    def start_stage(self, stage_name):
        pass

    def end_stage(self, stage_name, exit_status=0):
        pass

    def end_pipeline(self):
        pass

    def fail_pipeline(self, message='', **kwargs):
        pass

    @staticmethod
    def _format_error_message(message='', stacktrace=None):
        msg = 'Run failed.' + message
        if stacktrace:
            msg += '\nStack trace below:\n\n' + stacktrace
        return msg


default = NotificationCenter()
