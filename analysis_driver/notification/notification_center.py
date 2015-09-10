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

    def add_subscribers(self, run_id, **ntf_cfgs):
        for k in ntf_cfgs:
            notifier, config = ntf_cfgs[k]
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


default = NotificationCenter()
