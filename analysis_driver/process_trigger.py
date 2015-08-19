__author__ = 'mwham'
import argparse
import os
import logging
import logging.config

from analysis_driver import driver
from analysis_driver.config import default as cfg


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '--report',
        action='store_true',
        help='don\'t execute anything, report on status of datasets'
    )
    parser.add_argument(
        '--skip',
        help='mark a dataset as completed'
    )
    parser.add_argument(
        '--reset',
        help='unmark a dataset as complete/in progress for rerunning'
    )
    parser.add_argument(
        '--debug',
        action='store_true',
        help='override pipeline log level to debug'
    )

    args = parser.parse_args()

    logging.config.dictConfig(cfg.logging_config(debug=args.debug))
    logger = logging.getLogger(__name__)

    if args.report:
        report()
    elif args.skip:
        skip(args.skip)
    elif args.reset:
        reset(args.reset)
    else:
        datasets = scan_datasets()
        # Only process the first new dataset to be found. The rest will need to wait until the next scan
        if datasets:
            trigger(datasets[0])
        else:
            logger.info('No new datasets found')


def report():
    all_datasets = scan_datasets(ready_only=False)

    print('========= Process Trigger report =========')
    print('=== ready datasets ===')
    for d in all_datasets:
        if is_ready(d):
            print(d)

    print('=== unready datasets (RTA not complete) ===')
    for d in all_datasets:
        if is_unprocessed(d) and not rta_complete(d):
            print(d)

    print('=== active datasets ===')
    for d in all_datasets:
        if is_active(d):
            print(d)

    print('=== complete datasets ===')
    for d in all_datasets:
        if is_complete(d):
            print(d)


def scan_datasets(ready_only=True):
    all_datasets = [x for x in os.listdir(cfg['input_dir']) if os.path.isdir(os.path.join(cfg['input_dir'], x))]
    if ready_only:
        return [
            x for x in all_datasets
            if is_ready(x)
        ]
    else:
        return all_datasets


def trigger(dataset):
    active_lock = lock_file(dataset, 'active')
    complete_lock = lock_file(dataset, 'complete')
    assert not is_active(dataset)
    assert not is_complete(dataset)

    touch(active_lock)
    driver.pipeline(os.path.join(cfg['input_dir'], dataset))

    os.remove(active_lock)
    touch(complete_lock)


def skip(dataset):
    active_lock = lock_file(dataset, 'active')
    if os.path.isfile(active_lock):
        os.remove(active_lock)
    touch(lock_file(dataset, 'complete'))


def reset(dataset):
    for l in [
        lock_file(dataset, 'active'),
        lock_file(dataset, 'complete'),
        os.path.join(cfg['jobs_dir'], dataset, '.bcl2fastq_complete')
    ]:
        if os.path.isfile(l):
            os.remove(l)


def is_unprocessed(dataset):
    if not is_active(dataset) and not is_complete(dataset):
        return True
    else:
        return False


def is_ready(dataset):
    if is_unprocessed(dataset) and rta_complete(dataset):
        return True
    else:
        return False


def is_active(dataset):
    return os.path.isfile(lock_file(dataset, 'active'))


def is_complete(dataset):
    return os.path.isfile(lock_file(dataset, 'complete'))


def rta_complete(dataset):
    return os.path.isfile(os.path.join(cfg['input_dir'], dataset, 'RTAComplete.txt'))


def lock_file(dataset, status):
    return os.path.join(cfg['input_dir'], '.' + dataset + '.' + status)


def touch(file):
    open(file, 'w').close()
