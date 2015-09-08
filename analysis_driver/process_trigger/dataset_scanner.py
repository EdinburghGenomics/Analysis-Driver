__author__ = 'mwham'
import os
from analysis_driver.config import default as cfg


def report(all_datasets=False):
    datasets = scan_datasets()

    print('========= Process Trigger report =========')
    print('=== ready datasets ===')
    print('\n'.join(_fetch_by_status(datasets, 'ready')))

    print('=== unready datasets (RTA not complete) ===')
    print('\n'.join(_fetch_by_status(datasets, 'unready')))

    print('=== active datasets ===')
    print('\n'.join(_fetch_by_status(datasets, 'active')))

    print('=== complete datasets ===')
    complete_datasets = _fetch_by_status(datasets, 'complete')
    if all_datasets:
        pass
    elif complete_datasets == ['none']:
        pass
    else:
        complete_datasets = ['completed datasets are present', 'use --report-all to show']
    print('\n'.join(complete_datasets))


def scan_datasets():
    all_datasets = []
    for d in os.listdir(cfg['input_dir']):
        if os.path.isdir(os.path.join(cfg['input_dir'], d)):
            all_datasets.append((d, _status(d)))
            
    all_datasets.sort()
    return all_datasets


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


def _fetch_by_status(all_datasets, status):
    datasets = [d for (d, s) in all_datasets if s == status]
    if datasets:
        return datasets
    else:
        return ['none']


def _status(dataset):
    if is_ready(dataset):
        return 'ready'
    elif is_not_ready(dataset):
        return 'unready'
    elif _is_active(dataset):
        return 'active'
    elif _is_complete(dataset):
        return 'complete'
    else:
        return 'unknown'


def _is_processed(d):
    if _is_active(d) or _is_complete(d):
        return True
    else:
        return False


def is_ready(d):
    if not _is_processed(d) and _rta_complete(d):
        return True
    else:
        return False


def is_not_ready(d):
    if not _is_processed(d) and not _rta_complete(d):
        return True
    else:
        return False


def _is_active(d):
    return os.path.isfile(lock_file(d, 'active'))


def _is_complete(d):
    return os.path.isfile(lock_file(d, 'complete'))


def _rta_complete(d):
    return os.path.isfile(os.path.join(cfg['input_dir'], d, 'RTAComplete.txt'))


def lock_file(d, status):
    return os.path.join(cfg['input_dir'], '.' + d + '.' + status)


def touch(file):
    open(file, 'w').close()
