__author__ = 'mwham'
import argparse
import os
import logging
from sys import stdout
from analysis_driver.config import default as cfg
from analysis_driver.config import logging_default as log_cfg


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
    parser.add_argument(
        '--quiet',
        action='store_true',
        help='Do not log stdout'
    )

    args = parser.parse_args()

    if args.report:
        report()
    elif args.skip:
        skip(args.skip)
    elif args.reset:
        reset(args.reset)
    else:
        for d, s in scan_datasets():
            if s == 'ready':
                setup_run(d)
                setup_logging(d, args)
                logger = logging.getLogger('trigger')
                logger.info('Using config file at ' + cfg.config_file.name)
                logger.info('Triggering for dataset: ' + d)
                trigger(d)
                # only process the first new dataset found. The rest will need to wait until the next scan
                logger.info('Done')
                return 0

        setup_logging(None, args)
        logger = logging.getLogger('trigger')
        logger.debug('No new datasets found')

    return 0


def report():
    datasets = scan_datasets()

    print('========= Process Trigger report =========')
    print('=== ready datasets ===')
    print('\n'.join(_fetch_by_status(datasets, 'ready')))

    print('=== unready datasets (RTA not complete) ===')
    print('\n'.join(_fetch_by_status(datasets, 'unready')))

    print('=== active datasets ===')
    print('\n'.join(_fetch_by_status(datasets, 'active')))

    print('=== complete datasets ===')
    print('\n'.join(_fetch_by_status(datasets, 'complete')))


def scan_datasets():
    all_datasets = []
    for d in os.listdir(cfg['input_dir']):
        if os.path.isdir(os.path.join(cfg['input_dir'], d)):
            all_datasets.append((d, _status(d)))
            
    all_datasets.sort()
    return all_datasets


def trigger(dataset):
    active_lock = lock_file(dataset, 'active')
    complete_lock = lock_file(dataset, 'complete')
    assert not _is_active(dataset) and not _is_complete(dataset)

    touch(active_lock)
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


def setup_logging(dataset, args):
    if args.debug:
        log_cfg.log_level = logging.DEBUG
    else:
        log_cfg.log_level = logging.INFO

    # add stdout StreamHandler
    if not args.quiet:
        log_cfg.add_handler('stdout', logging.StreamHandler(stream=stdout))  # TODO: do tthis via the config instead

    # add dataset-specific FileHandler
    if dataset:
        log_cfg.add_handler(
            'file',
            logging.FileHandler(
                filename=os.path.join(
                    cfg['jobs_dir'],
                    dataset,
                    'analysis_driver.log'
                )
            )
        )

    # add user-defined handlers
    if cfg['logging']['handlers']:
        for name, info in cfg['logging']['handlers'].items():
            log_cfg.add_handler(name, logging.FileHandler(info['filename']))

    log_cfg.lock()


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
    if _is_ready(dataset):
        return 'ready'
    elif not _is_processed(dataset) and not _rta_complete(dataset):
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


def _is_ready(d):
    if not _is_processed(d) and _rta_complete(d):
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
