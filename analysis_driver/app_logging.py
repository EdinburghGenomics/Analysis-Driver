import logging
import logging.config
import logging.handlers
from analysis_driver.config import default as cfg


class LoggingConfiguration:
    """Stores Loggers, Formatters and Handlers"""
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

    def get_logger(self, name, level=None):
        """
        Return a logging.Logger object with formatters and handlers added.
        :param name: A name to assign to the logger (usually __name__)
        :rtype: logging.Logger
        """
        logger = logging.getLogger(name)
        if level is None:
            level = self.log_level
        logger.setLevel(level)
        for h in self.handlers:
            logger.addHandler(h)
        self.loggers.add(logger)
        return logger

    def add_handler(self, handler, level=logging.NOTSET):
        """Add a handler, set its format/level if needed and register all loggers to it"""
        if level == logging.NOTSET:
            level = self.log_level
        handler.setLevel(level)
        handler.setFormatter(self.formatter)
        for l in self.loggers:
            l.addHandler(handler)
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

    def configure_handlers_from_config(self, handler_cfg):
        if not handler_cfg:
            return None

        configurator = logging.config.BaseConfigurator({})
        handler_classes = {
            'stream_handlers': logging.StreamHandler,
            'file_handlers': logging.FileHandler,
            'timed_rotating_file_handlers': logging.handlers.TimedRotatingFileHandler
        }

        for handler_type in handler_classes:
            for handler_cfg in handler_cfg.get(handler_type, []):
                level = logging.getLevelName(handler_cfg.pop('level', self.log_level))

                if 'stream' in handler_cfg:
                    handler_cfg['stream'] = configurator.convert(handler_cfg['stream'])
                handler = handler_classes[handler_type](**handler_cfg)
                self.add_handler(handler, level)


logging_default = LoggingConfiguration()


class AppLogger:
    """
    Mixin class for logging. An object subclassing this can log using its class name. Contains a
    logging.Logger object and exposes its log methods.
    """
    log_cfg = logging_default
    __logger = None

    def debug(self, msg, *args):
        self._logger.debug(msg, *args)

    def info(self, msg, *args):
        self._logger.info(msg, *args)

    def warning(self, msg, *args):
        self._logger.warning(msg, *args)

    def error(self, msg, *args):
        self._logger.error(msg, *args)

    def critical(self, msg, *args):
        self._logger.critical(msg, *args)

    @property
    def _logger(self):
        if self.__logger is None:
            self.__logger = self.log_cfg.get_logger(self.__class__.__name__)
        return self.__logger
