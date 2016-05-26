from .executor import Executor, StreamExecutor, ClusterExecutor, ArrayExecutor
from analysis_driver.config import default as cfg
from analysis_driver.app_logging import logging_default as log_cfg

app_logger = log_cfg.get_logger('executor')


def execute(cmds, env=None, **kwargs):
    """
    :param list[str] cmds: A list where each item is a list of strings to be passed to Executor
    :param bool cluster:
    :param kwargs:
    :return: Executor
    """
    if env is None:
        env = cfg.get('job_execution', 'local')

    if env == 'local':
        e = ArrayExecutor(cmds, stream=kwargs.get('stream', False))
    else:
        e = ClusterExecutor(cmds, qsub=cfg.query('tools', 'qsub', ret_default='qsub'), **kwargs)

    e.start()
    return e
