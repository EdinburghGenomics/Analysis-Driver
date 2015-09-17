__author__ = 'mwham'
import os.path
from time import sleep
from analysis_driver.dataset_scanner import dataset_status
from analysis_driver.executor import StreamExecutor
from analysis_driver.app_logging import get_logger

app_logger = get_logger('tt_agent')


def transfer_to_int_dir(dataset, from_dir, to_dir, repeat_delay):
    """
    rsync -aqu --size-only --partial from_dir/dataset to_dir
    :param str dataset:
    :param str from_dir:
    :param str to_dir:
    :param int repeat_delay:
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
