from .executor import SimpleExecutor, StreamExecutor, ClusterExecutor, ArrayExecutor
from analysis_driver.config import default as cfg
from analysis_driver.app_logging import get_logger

app_logger = get_logger('executor')


def execute(cmds, env=None, **kwargs):
    """
    :param list[str] cmds: A list where each item is a list of strings to be passed to SimpleExecutor
    :param bool cluster:
    :param kwargs:
    :return: Executor
    """
    if env is None:
        env = cfg.get('job_execution', 'local')

    if env == 'local':
        e = ArrayExecutor(cmds, stream=kwargs.get('stream', False))
    else:
        e = ClusterExecutor(cmds, **kwargs)

    e.start()
    return e
