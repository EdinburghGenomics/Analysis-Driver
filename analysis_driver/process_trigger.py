__author__ = 'mwham'
import os
from time import sleep
from analysis_driver.dataset_scanner import DATASET_READY, DATASET_NEW
from analysis_driver import executor
from analysis_driver.app_logging import get_logger
from analysis_driver.config import default as cfg


app_logger = get_logger('proctrigger')


def trigger(dataset):
    """
    Decide whether to rsync a dataset to an intermediate dir and run driver.pipeline on it.
    :param str dataset: A dataset id
    """
    status = dataset.dataset_status

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

    dataset.start()

    from analysis_driver import driver
    exit_status = driver.pipeline(dataset_dir, dataset)

    if exit_status != 0:
        dataset.fail()
    else:
        dataset.succeed()

    return exit_status


def _transfer_to_int_dir(dataset, from_dir, to_dir, repeat_delay):
    """
    rsync -aqu --size-only --partial from_dir/dataset to_dir
    """
    exit_status = 0
    app_logger.info('Starting transfer')
    rsync_cmd = 'rsync -aqu --size-only --partial %s %s' % (os.path.join(from_dir, dataset.name), to_dir)

    while dataset.dataset_status != DATASET_READY:
        exit_status += executor.execute([rsync_cmd], job_name='rsync', run_id=dataset.name, walltime=36).join()
        sleep(repeat_delay)



    # one more rsync after the RTAComplete is created. After this, everything should be synced
    sleep(repeat_delay)
    exit_status += executor.execute([rsync_cmd], job_name='rsync', run_id=dataset.name, walltime=36).join()
    assert os.path.isfile(os.path.join(dataset.path, 'RTAComplete.txt'))
    app_logger.info('Transfer complete with exit status ' + str(exit_status))
