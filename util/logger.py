__author__ = 'mwham'

import logging


class AppLogger:
    __logger = None

    def log(self, msg, level='INFO', error_class=None):
        if self.__logger is None:
            self.__logger = logging.getLogger(self.__class__.__name__)

        if level == 'DEBUG':
            self.__logger.debug(msg)
        elif level == 'INFO':
            self.__logger.info(msg)
        elif level == 'WARN':
            self.__logger.warn(msg)
        elif level == 'ERROR':
            self.__logger.error(msg)
        elif level == 'CRITICAL':
            if error_class:
                raise error_class(msg)
            else:
                self.__logger.critical(msg)
        else:
            raise ValueError('Invalid log level: ' + level)
