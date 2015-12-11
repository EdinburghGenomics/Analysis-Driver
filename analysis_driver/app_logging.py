__author__ = 'mwham'
import logging
from analysis_driver.config import logging_default as log_cfg


class AppLogger:
    """
    Mixin class for logging. An object subclassing this can log directly from itself. Contains a
    logging.Logger object, which is controlled by public methods. All methods are voids, i.e. return None.
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
        self._logger.warning(msg)

    def error(self, msg):
        self._check_logger()
        self._logger.error(msg)

    def critical(self, msg, error_class=None):
        """
        Log at the level logging.CRITICAL. Can also raise an error if an error_class is given.
        :param error_class: An Error to be raised, e.g. ValueError, AssertionError
        """
        self._check_logger()

        if error_class:
            raise error_class(msg)
        else:
            self._logger.critical(msg)

    def _check_logger(self):
        if self._logger is None:
            self._logger = get_logger(self.__class__.__name__)


def get_logger(name):
    """
    Return a logging.Logger object with formatters and handlers added.
    :param name: A name to assign to the logger (usually __name__)
    :rtype: logging.Logger
    """
    logger = logging.getLogger(name)
    logger.setLevel(log_cfg.default_level)
    for name, h in log_cfg.handlers.items():
        logger.addHandler(h)

    return logger
