__author__ = 'mwham'
import os
from analysis_driver.config import default as cfg
from .dataset_scanner import report, skip, reset, scan_datasets, lock_file, touch
from . import tt_agent


def trigger(dataset):
    active_lock = lock_file(dataset, 'active')
    complete_lock = lock_file(dataset, 'complete')

    touch(active_lock)
    if cfg.get('raw_dir'):
        tt_agent.rsync(dataset)

    from analysis_driver import driver
    driver.pipeline(os.path.join(cfg['input_dir'], dataset))

    os.remove(active_lock)
    touch(complete_lock)


def setup_run(dataset):
    for d in ['fastq_dir', 'jobs_dir']:
        try:
            os.makedirs(os.path.join(cfg[d], dataset))
        except FileExistsError:
            pass
