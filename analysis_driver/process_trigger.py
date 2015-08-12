__author__ = 'mwham'
import argparse
import os
import logging
import logging.config

from analysis_driver import driver
from analysis_driver.config import default as cfg
from analysis_driver.util import ProcessTriggerError


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '--report',
        action='store_true',
        help='Don\'t execute anything, report on status of datasets'
    )
    parser.add_argument(
        '--debug',
        action='store_true',
        help='Override pipeline log level to debug'
    )

    args = parser.parse_args()

    logging.config.dictConfig(cfg.logging_config(debug=args.debug))
    logger = logging.getLogger(__name__)

    pt_lock_file = os.path.join(cfg['input_dir'], '.proctrigger.lock', 'w')
    if os.path.isfile(pt_lock_file):
        raise ProcessTriggerError(
            'Lock file present. Is another Process Trigger running? Check for ' + pt_lock_file
        )
    if args.report:
        report()
    else:
        datasets = scan_datasets()
        # Only process the first new dataset to be found. The rest will need to wait until the next scan
        if datasets:
            trigger(datasets[0])
        else:
            logger.info('No new datasets found')


def report():
    print('========= Process Trigger report =========\n')
    print('=== ready datasets ===')
    for dataset in scan_datasets(get_all=True):
        if is_ready(dataset):
            print(dataset)

    print('=== unready datasets (RTA not complete) ===')
    for dataset in scan_datasets(get_all=True):
        if is_unprocessed(dataset) and not rta_complete(dataset):
            print(dataset)

    print('=== active datasets ===')
    for dataset in scan_datasets(get_all=True):
        if is_active(dataset):
            print(dataset)

    print('=== complete datasets ===')
    for dataset in scan_datasets(get_all=True):
        if is_complete(dataset):
            print(dataset)


def scan_datasets(get_all=False):
    all_datasets = [x for x in os.listdir(cfg['input_dir']) if os.path.isdir(os.path.join(cfg['input_dir'], x))]
    if get_all:
        return all_datasets
    else:
        return [
            x for x in all_datasets
            if is_ready(x)
        ]


def is_unprocessed(dataset):
    if not is_active(dataset) and not is_complete(dataset):
        return True
    else:
        return False


def is_active(dataset):
    return os.path.isfile(lock_file(dataset, 'active'))


def is_complete(dataset):
    return os.path.isfile(lock_file(dataset, 'complete'))


def is_ready(dataset):
    if is_unprocessed(dataset) and rta_complete(dataset):
        return True
    else:
        return False


def lock_file(dataset, status):
    return os.path.join(cfg['input_dir'], '.' + dataset + '.' + status)


def rta_complete(dataset):
    return os.path.isfile(os.path.join(cfg['input_dir'], dataset, 'RTAComplete.txt'))


def trigger(dataset):
    active_lock = lock_file(dataset, 'active')
    complete_lock = lock_file(dataset, 'complete')
    assert not is_active(dataset)
    assert not is_complete(dataset)

    touch(active_lock)
    driver.pipeline(os.path.join(cfg['input_dir'], dataset))

    os.remove(active_lock)
    touch(complete_lock)


def touch(file):
    open(file, 'w').close()
