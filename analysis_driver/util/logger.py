__author__ = 'mwham'
import logging


class AppLogger:
    """
    Mixin class for logging. Contains a logging.Logger object, which is controlled by public methods. All
    methods are voids, i.e. return None.

    Public methods:
        All of these tell self._logger to log at the corresponding level.
        debug()
        info()
        warn()
        error()
        critical()
    """
    _logger = None

    def debug(self, msg):
        self._check_logger()
        self._logger.debug(msg)

    def info(self, msg):
        self._check_logger()
        self._logger.info(msg)

    def warn(self, msg):
        self._check_logger()
        self._logger.warn(msg)

    def error(self, msg):
        self._check_logger()
        self._logger.error(msg)

    def critical(self, msg, error_class=None):
        """
        Log at the level logging.CRITICAL. Can also raise an error if an error_class is given.
        :param msg: Log message
        :param error_class: An Error to be raised, e.g. ValueError, AssertionError
        :return: None
        :raises: An arbitrary Error class if error_class is specified
        """
        self._check_logger()

        if error_class:
            raise error_class(msg)
        else:
            self._logger.critical(msg)

    def _check_logger(self):
        if self._logger is None:
            self._logger = logging.getLogger(self.__class__.__name__)


class NamedAppLogger(AppLogger):
    """
    Subclass of AppLogger, which can be created with a name. Useful for logging from non-object-oriented
    modules.
    """
    def __init__(self, name):
        self._logger = logging.getLogger(name)
