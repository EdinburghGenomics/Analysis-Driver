__author__ = 'mwham'
import os
from analysis_driver.config import default as cfg


def report(all_datasets=False):
    datasets = scan_datasets()

    print('========= Process Trigger report =========')
    print('=== new datasets ===')
    print('\n'.join(_fetch_by_status(datasets, 'new', 'new, rta complete')))
    print('=== transferring datasets ===')
    print('\n'.join(_fetch_by_status(datasets, 'transferring', 'transferring, rta complete')))

    print('=== active datasets ===')
    print('\n'.join(_fetch_by_status(datasets, 'active')))

    if all_datasets:
        print('=== complete datasets ===')
        print('\n'.join(_fetch_by_status(datasets, 'complete')))
        print('=== failed datasets ===')
        print('\n'.join(_fetch_by_status(datasets, 'failed')))
        print('=== aborted datasets ===')
        print('\n'.join(_fetch_by_status(datasets, 'aborted')))
    else:
        print('=== other datasets ===')
        print('\n'.join(['completed datasets are present', 'use --report-all to show']))


def scan_datasets():
    triggerignore_file = os.path.join(
        cfg.get('lock_file_dir', cfg['input_dir']),
        '.triggerignore'
    )
    if os.path.isfile(triggerignore_file):
        triggerignore = [x.strip() for x in open(triggerignore_file).readlines()]
    else:
        triggerignore = []

    all_datasets = []
    for d in os.listdir(cfg['input_dir']):
        if os.path.isdir(os.path.join(cfg['input_dir'], d)) and d not in triggerignore:
            all_datasets.append((d, dataset_status(d)))
            
    all_datasets.sort()
    return all_datasets


def skip(dataset):
    reset(dataset)
    touch(lock_file(dataset, 'complete'))


def reset(dataset):
    for s in ('active', 'complete', 'transferring', 'aborted', 'failed'):
        _rm(lock_file(dataset, s))


def abort(dataset):
    reset(dataset)
    touch(lock_file(dataset, 'aborted'))


def switch_status(dataset, status):
    reset(dataset)
    touch(lock_file(dataset, status))


def _fetch_by_status(all_datasets, *statuses):
    datasets = [d for (d, s) in all_datasets if s in statuses]
    if datasets:
        return datasets
    else:
        return ['none']
# TODO: remove duplicate code!

def dataset_status(dataset):
    transferring = _transferring(dataset)
    rta_complete = _rta_complete(dataset)

    if _complete(dataset):
        assert rta_complete
        return 'complete'
    elif _active(dataset):
        assert rta_complete
        return 'active'
    elif _aborted(dataset):
        return 'aborted'
    elif _failed(dataset):
        return 'failed'

    elif transferring and rta_complete:
        return 'transferring, rta complete'
    elif transferring:
        return 'transferring'
    elif rta_complete:
        return 'new, rta complete'
    
    else:
        return 'new'


def _active(dataset):
    if os.path.isfile(lock_file(dataset, 'active')):
        for status in (_aborted, _complete, _transferring, _failed):
            assert not status(dataset)
        return True
    else:
        return False


def _complete(dataset):
    if os.path.isfile(lock_file(dataset, 'complete')):
        for status in (_active, _aborted, _transferring, _failed):
            assert not status(dataset)
        return True
    else:
        return False


def _transferring(dataset):
    if os.path.isfile(lock_file(dataset, 'transferring')):
        for status in (_active, _complete, _aborted, _failed):
            assert not status(dataset)
        return True
    else:
        return False


def _aborted(dataset):
    if os.path.isfile(lock_file(dataset, 'aborted')):
        for status in (_active, _complete, _transferring, _failed):
            assert not status(dataset)
        return True
    else:
        return False


def _failed(dataset):
    if os.path.isfile(lock_file(dataset, 'failed')):
        for status in (_active, _complete, _transferring, _aborted):
            assert not status(dataset)
        return True
    else:
        return False


def _rta_complete(dataset):
    return os.path.isfile(os.path.join(cfg['input_dir'], dataset, 'RTAComplete.txt'))


def lock_file(dataset, status):
    return os.path.join(cfg.get('lock_file_dir', cfg['input_dir']), '.' + dataset + '.' + status)


def touch(file):
    open(file, 'w').close()


def _rm(file):
    if os.path.isfile(file):
        os.remove(file)
