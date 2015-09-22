__author__ = 'mwham'
import os
from time import sleep
from analysis_driver.dataset_scanner import dataset_status, switch_status
from analysis_driver.executor import StreamExecutor
from analysis_driver.app_logging import get_logger
from analysis_driver.config import default as cfg


app_logger = get_logger('proctrigger')


def trigger(dataset, use_intermediate_dir):
    """
    Decide whether to rsync a dataset to an intermediate dir and run driver.pipeline on it.
    :param str dataset: A dataset id
    :param use_intermediate_dir: Whether to rsync to an intermediate dir
    """
    if use_intermediate_dir:
        assert dataset_status(dataset) in ('new', 'new, rta complete')
        switch_status(dataset, 'transferring')
        transfer_to_int_dir(
            dataset,
            cfg['input_dir'],
            cfg['intermediate_dir'],
            repeat_delay=int(cfg['tt_agent_delay'])
        )
        dataset_dir = cfg['intermediate_dir']
    else:
        assert dataset_status(dataset) == 'new, rta complete'
        dataset_dir = cfg['input_dir']

    switch_status(dataset, 'active')
    from analysis_driver import driver
    exit_status = driver.pipeline(os.path.join(dataset_dir, dataset))

    if exit_status != 0:
        switch_status(dataset, 'failed')
    else:
        switch_status(dataset, 'complete')

    return exit_status


def transfer_to_int_dir(dataset, from_dir, to_dir, repeat_delay):
    """
    rsync -aqu --size-only --partial from_dir/dataset to_dir
    """
    app_logger.info('Starting transfer')
    while dataset_status(dataset) != 'transferring, rta complete':
        _run(['rsync', '-aqu', '--size-only', '--partial', os.path.join(from_dir, dataset), to_dir])
        sleep(repeat_delay)

    # one more rsync after the RTAComplete is created. After this, everything should be synced
    sleep(repeat_delay)
    _run(['rsync', '-aqu', '--size-only', '--partial', os.path.join(from_dir, dataset), to_dir])
    assert os.path.isfile(os.path.join(to_dir, dataset, 'RTAComplete.txt'))
    app_logger.info('Transfer complete')


def _run(cmd):
    executor = StreamExecutor(cmd)
    executor.run()
