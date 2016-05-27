from .executor import Executor, StreamExecutor, ArrayExecutor, PBSExecutor, SlurmExecutor
from analysis_driver.config import default as cfg
from analysis_driver.exceptions import AnalysisDriverError


def local_execute(*cmds, parallel=True):
    if len(cmds) == 1:
        if parallel:
            e = StreamExecutor(cmds[0])
        else:
            e = Executor(cmds[0])
    else:
        e = ArrayExecutor(cmds, stream=parallel)

    e.start()
    return e


def cluster_execute(*cmds, env=None, prelim_cmds=None, **cluster_config):
    if env == 'pbs':
        cls = PBSExecutor
    elif env == 'slurm':
        cls = SlurmExecutor
    else:
        raise AnalysisDriverError('Unknown execution environment: ' + env)

    e = cls(*cmds, prelim_cmds=prelim_cmds, **cluster_config)
    e.start()
    return e


def execute(*cmds, env=None, prelim_cmds=None, **cluster_config):
    """
    :param list[str] cmds: A list where each item is a list of strings to be passed to Executor
    :param str env:
    :param cluster_config:
    :return: Executor
    """
    if env is None:
        env = cfg.get('job_execution', 'local')

    if env == 'local':
        return local_execute(*cmds)
    else:
        return cluster_execute(*cmds, env=env, prelim_cmds=prelim_cmds, **cluster_config)
