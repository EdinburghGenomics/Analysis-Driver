__author__ = 'mwham'
import os
from time import sleep
from analysis_driver.dataset_scanner import RunScanner, DATASET_READY, DATASET_NEW
from analysis_driver import executor
from analysis_driver.app_logging import get_logger
from analysis_driver.config import default as cfg


app_logger = get_logger('proctrigger')

scanner = RunScanner(cfg)

def trigger(dataset):
    """
    Decide whether to rsync a dataset to an intermediate dir and run driver.pipeline on it.
    :param str dataset: A dataset id
    """
    status = scanner.dataset_status(dataset)

    if cfg.get('intermediate_dir'):
        assert status in [DATASET_NEW, DATASET_READY], 'Invalid dataset status: ' + status
        _transfer_to_int_dir(
            dataset,
            cfg['input_dir'],
            cfg['intermediate_dir'],
            repeat_delay=int(cfg.get('tt_agent_delay', 120))
        )
        dataset_dir = cfg['intermediate_dir']
    else:
        assert status in [DATASET_READY], 'Invalid dataset status: ' + status
        dataset_dir = cfg['input_dir']

    scanner.start(dataset)
    from analysis_driver import driver
    exit_status = driver.pipeline(os.path.join(dataset_dir, dataset))

    if exit_status != 0:
        scanner.fail(dataset)
    else:
        scanner.succeed(dataset)

    return exit_status


def _transfer_to_int_dir(dataset, from_dir, to_dir, repeat_delay):
    """
    rsync -aqu --size-only --partial from_dir/dataset to_dir
    """
    exit_status = 0
    app_logger.info('Starting transfer')
    rsync_cmd = 'rsync -aqu --size-only --partial %s %s' % (os.path.join(from_dir, dataset), to_dir)

    while scanner.dataset_status(dataset) != DATASET_READY:
        exit_status += executor.execute([rsync_cmd], job_name='rsync', run_id=dataset, walltime=36).join()
        sleep(repeat_delay)

    # one more rsync after the RTAComplete is created. After this, everything should be synced
    sleep(repeat_delay)
    exit_status += executor.execute([rsync_cmd], job_name='rsync', run_id=dataset, walltime=36).join()
    assert os.path.isfile(os.path.join(to_dir, dataset, 'RTAComplete.txt'))
    app_logger.info('Transfer complete with exit status ' + str(exit_status))
