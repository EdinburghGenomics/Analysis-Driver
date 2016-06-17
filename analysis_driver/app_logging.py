from egcg_core.app_logging import LoggingConfiguration, AppLogger as EGCGLogger
from analysis_driver.config import default as cfg


log_cfg = LoggingConfiguration(cfg)


class AppLogger(EGCGLogger):
    log_cfg = log_cfg
