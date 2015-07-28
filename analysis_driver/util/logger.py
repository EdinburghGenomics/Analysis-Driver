__author__ = 'mwham'
import logging


class AppLogger:
    __logger = None

    def debug(self, msg):
        self._check_logger()
        self.__logger.debug(msg)

    def info(self, msg):
        self._check_logger()
        self.__logger.info(msg)

    def warn(self, msg):
        self._check_logger()
        self.__logger.warn(msg)

    def error(self, msg):
        self._check_logger()
        self.__logger.error(msg)

    def critical(self, msg, error_class=None):
        self._check_logger()

        if error_class:
            raise error_class(msg)
        else:
            self.__logger.critical(msg)

    def _check_logger(self):
        if self.__logger is None:
            self.__logger = logging.getLogger(self.__class__.__name__)
