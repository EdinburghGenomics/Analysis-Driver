__author__ = 'mwham'
import os
import argparse
from config import default as cfg
from process_trigger import ProcessTriggerError


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '--report',
        action='store_true',
        help='Don\'t execute anything, report on status of datasets'
    )

    args = parser.parse_args()

    pt_lock_file = os.path.join(cfg['input_dir'], '.proctrigger.lock', 'w')
    if os.path.isfile(pt_lock_file):
        raise ProcessTriggerError(
            'Lock file present. Is another Process Trigger running? Check for ' + pt_lock_file
        )
    else:
        open(pt_lock_file, 'w').close()

    if args.report:
        report()
    else:
        for dataset in scan_datasets():
            pass
            # set off a new analysis_driver.main for each dataset


def report():
    print('========= Process Trigger report =========\n')
    print('=== ready datasets ===')
    for dataset in os.listdir(cfg['input_dir']):
        if is_ready(dataset):
            print(dataset)

    print('=== unready datasets (RTA not complete) ===')
    for dataset in os.listdir(cfg['input_dir']):
        if is_unprocessed(dataset) and not rta_complete(dataset):
            print(dataset)

    print('=== active datasets ===')
    for dataset in os.listdir(cfg['input_dir']):
        if is_active(dataset):
            print(dataset)

    print('=== complete datasets ===')
    for dataset in os.listdir(cfg['input_dir']):
        if is_complete(dataset):
            print(dataset)


def scan_datasets():
    return [os.path.join(cfg['input_dir'], x) for x in os.listdir(cfg['input_dir']) if is_ready(x)]


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
    pass
    # touch an 'active' lock file, set off an analysis_driver, then touch a 'complete' lock file
    # need to move the executor module to the top level
