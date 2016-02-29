__author__ = 'mwham'
import logging
import logging.config
import logging.handlers
from analysis_driver.config import default as cfg


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
    logger.setLevel(log_cfg.log_level)
    for h in log_cfg.handlers:
        logger.addHandler(h)
    return logger


class LoggingConfiguration:
    """Stores logging Formatters and Handlers."""
    def __init__(self):
        self.default_formatter = logging.Formatter(
            fmt=cfg.query('logging', 'format', ret_default='[%(asctime)s][%(name)s][%(levelname)s] %(message)s'),
            datefmt=cfg.query('logging', 'datefmt', ret_default='%Y-%b-%d %H:%M:%S')
        )
        self.blank_formatter = logging.Formatter()
        self.formatter = self.default_formatter
        self.handlers = set()
        self.loggers = set()
        self.log_level = logging.INFO

    def add_handler(self, handler):
        """
        :param logging.FileHandler handler:
        """
        if handler.level == logging.NOTSET:
            handler.setLevel(self.log_level)
        handler.setFormatter(self.formatter)
        self.handlers.add(handler)

    def set_log_level(self, level):
        self.log_level = level
        for h in self.handlers:
            h.setLevel(self.log_level)
        for l in self.loggers:
            l.setLevel(self.log_level)

    def set_formatter(self, formatter):
        """Set all handlers to use formatter"""
        self.formatter = formatter
        for h in self.handlers:
            h.setFormatter(self.formatter)

    def configure_handlers_from_config(self, cfg_logging):
        configurator = logging.config.BaseConfigurator({})
        handler_classes = {
            'stream_handlers': logging.StreamHandler,
            'file_handlers': logging.FileHandler,
            'timed_rotating_file_handlers': logging.handlers.TimedRotatingFileHandler
        }

        for handler_type in handler_classes:
            for handler_cfg in cfg_logging.get(handler_type, []):
                level = logging.getLevelName(handler_cfg.pop('level', 'WARNING'))
                if 'stream' in handler_cfg:
                    handler_cfg['stream'] = configurator.convert(handler_cfg['stream'])
                handler = handler_classes[handler_type](**handler_cfg)
                handler.setLevel(level)
                self.add_handler(handler)


log_cfg = LoggingConfiguration()
