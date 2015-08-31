from analysis_driver.notification.email_notification import EmailNotification
from analysis_driver.config import default as cfg
from analysis_driver.notification.log_notification import LogNotification

__author__ = 'tcezard'

class _Method:
    # some magic to bind an notification center method to its subscribers.
    def __init__(self, send, name):
        self.__send = send
        self.__name = name
    def __call__(self, *args, **kwargs):
        return self.__send(self.__name, *args, **kwargs)


class NotificationCenter:
    '''Object dispatching notification to subscribers.
    The suscribers are defined in a config file'''
    def __init__(self, config):
        self.subscribers=[]
        if config:
            self._init_with_config(config['notification'])

    def _init_with_config(self,config):
        if 'email_notification' in config:
            self.subscribers.append(EmailNotification(config['email_notification']))
        if 'log_notification' in config:
            self.subscribers.append(LogNotification(config['log_notification']))

    def _pass_to_subscriber(self, function_name, *args, **kwargs):
        for subscriber in self.subscribers:
            f = getattr(subscriber, function_name)
            if f and callable(f):
                f(*args, **kwargs)
            else:
                raise(ValueError("Method {} does not exist in {}".format(function_name, subscriber.__name__)))

    def __getattr__(self, name):
        # magic method dispatcher
        return _Method(self._pass_to_subscriber, name)


default_notification_center = NotificationCenter(cfg)