__author__ = 'mwham'
import os
from glob import glob
from collections import defaultdict
from analysis_driver.app_logging import get_logger
from analysis_driver.config import default as cfg

app_logger = get_logger('scanner')


def report(all_datasets=False):
    datasets = scan_datasets()

    print('========= Process Trigger report =========')
    print('dataset location: ' + cfg['input_dir'])
    for status in ('new', 'new, rta complete', 'transferring', 'transferring, rta complete', 'active'):
        ds = datasets.pop(status, [])
        if ds:
            print('=== ' + status + ' ===')
            print('\n'.join((os.path.basename(d) for d in ds)))

    if any((datasets[s] for s in datasets)):
        if all_datasets:
            for status in sorted(datasets):
                print('=== ' + status + ' ===')
                print('\n'.join((os.path.basename(x) for x in datasets[status])))
        else:
            print('=== other datasets ===')
            print('\n'.join(('other datasets present', 'use --report-all to show')))

    print('_' * 42)


def scan_datasets():
    lock_dir = cfg.get('lock_file_dir', cfg['input_dir'])
    input_dir = cfg['input_dir']
    triggerignore = os.path.join(lock_dir, '.triggerignore')

    ignorables = []
    if os.path.isfile(triggerignore):
        with open(triggerignore, 'r') as f:
            for p in f.readlines():
                if not p.startswith('#'):
                    ignorables.extend(glob(os.path.join(input_dir, p.rstrip('\n'))))
    app_logger.debug('Ignoring %s datasets' % len(ignorables))

    n_datasets = 0
    datasets = defaultdict(list)
    for d in glob(os.path.join(input_dir, '*')):
        if os.path.isdir(d) and d not in ignorables:
            datasets[dataset_status(d)].append(d)
            n_datasets += 1

    app_logger.debug('Found %s datasets' % n_datasets)
    return datasets


def reset(dataset):
    _rm(*glob(lock_file(dataset, '*')))


def switch_status(dataset, status):
    reset(dataset)
    _touch(lock_file(dataset, status))


def dataset_status(dataset):
    dataset_lock_files = glob(lock_file(dataset, '*'))
    assert len(dataset_lock_files) < 2
    if dataset_lock_files:
        lf_status = dataset_lock_files[0].split('.')[-1]
    else:
        lf_status = 'new'

    rta_complete = _rta_complete(dataset)

    if lf_status in ('complete', 'active'):
        assert rta_complete
        return lf_status

    elif lf_status in ('aborted', 'failed'):
        return lf_status

    else:
        if rta_complete:
            lf_status += ', rta complete'
        return lf_status


def _rta_complete(dataset):
    return os.path.isfile(os.path.join(cfg['input_dir'], dataset, 'RTAComplete.txt'))


def lock_file(dataset, status):
    return os.path.join(
        cfg.get('lock_file_dir', cfg['input_dir']),
        '.' + os.path.basename(dataset) + '.' + status
    )


def _touch(file):
    open(file, 'w').close()


def _rm(*files):
    for f in files:
        if os.path.isfile(f):
            os.remove(f)
