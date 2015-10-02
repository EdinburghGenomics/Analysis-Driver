from .executor import Executor, StreamExecutor, ClusterExecutor
from analysis_driver import writer
from analysis_driver.app_logging import get_logger

app_logger = get_logger('executor')


def execute(cmds, cluster=False, **kwargs):
    """

    :param list cmds: A list where each item is a list of strings to be passed to Executor
    :param bool cluster:
    :param kwargs:
    :return:
    """
    if not cluster:
        return _local_execute(cmds, stream=kwargs.get('stream', False))
    else:
        return _cluster_execute(cmds, **kwargs)


def _local_execute(cmds, stream):

    for cmd in cmds:
        if stream:
            e = StreamExecutor(cmd)
            e.start()
            exit_status = e.join()
            app_logger.info(exit_status)
            return exit_status

        else:
            e = Executor(cmd)
            out, err = e.run()
            app_logger.info('Output:')
            for line in out.decode('utf-8').split('\n'):
                app_logger.info(line)
            for line in err.decode('utf-8').split('\n'):
                app_logger.error(line)

            return e.proc.poll()  # exit status


def _cluster_execute(cmds, **kwargs):
    prelim_cmds = kwargs.pop('prelim_cmds', None)
    w = writer.get_script_writer(jobs=len(cmds), **kwargs)
    w.write_jobs(cmds, prelim_cmds)

    e = ClusterExecutor(w.script_name)
    e.start()
    return e
