__author__ = 'mwham'
import argparse
import logging
import os
import sys
from analysis_driver.app_logging import get_logger
from analysis_driver.config import default as cfg, logging_default as log_cfg
from analysis_driver.notification import default as ntf, LogNotification, EmailNotification
from analysis_driver.exceptions import AnalysisDriverError
from analysis_driver.dataset_scanner import SampleScanner, DATASET_READY, DATASET_FORCE_READY, RunScanner


def main():
    args = _parse_args()

    if args.debug:
        log_cfg.default_level = logging.DEBUG

    logging_handlers = cfg.query('logging', 'handlers')
    if logging_handlers:
        for name, config in logging_handlers.items():
            if 'filename' in config:
                handler = logging.FileHandler(filename=config['filename'], mode='a')
            elif 'stream' in config:
                s = None
                if config['stream'] == 'ext://sys.stdout':
                    s = sys.stdout
                elif config['stream'] == 'ext://sys.stderr':
                    s = sys.stderr
                handler = logging.StreamHandler(stream=s)
            else:
                raise AnalysisDriverError('Invalid logging configuration: %s %s' % name, str(config))
            log_cfg.add_handler(name, handler, config.get('level', log_cfg.default_level))

    if args.run:
        if 'run' in cfg:
            cfg.merge(cfg['run'])
        scanner = RunScanner(cfg)
    else:
        assert args.sample
        if 'sample' in cfg:
            cfg.merge(cfg['sample'])
        scanner = SampleScanner(cfg)

    if any([args.abort, args.skip, args.reset, args.force, args.report, args.report_all]):
        for d in args.abort:
            scanner.get(d).abort()
        for d in args.skip:
            dataset = scanner.get(d)
            dataset.reset()
            dataset.start()
            dataset.succeed()
        for d in args.reset:
            scanner.get(d).reset()
        for d in args.force:
            scanner.get(d).force()

        if args.report:
            scanner.report()
        elif args.report_all:
            scanner.report(all_datasets=True)
        return 0

    all_datasets = scanner.scan_datasets()
    ready_datasets = all_datasets.get(DATASET_READY, []) + all_datasets.get(DATASET_FORCE_READY, [])

    if not ready_datasets:
        return 0
    else:
        # Only process the first new dataset found. Run through Cron, this will result in one new pipeline
        # being kicked off per minute.
        return _process_dataset(ready_datasets[0])


def setup_logging(d):
    log_repo = cfg.query('logging', 'repo')
    if log_repo:
        handler = logging.FileHandler(filename=os.path.join(log_repo, d.name + '.log'), mode='w')
        log_cfg.add_handler(d, handler, log_cfg.default_level)

    log_cfg.add_handler(
        'dataset',
        logging.FileHandler(filename=os.path.join(cfg['jobs_dir'], d.name, 'analysis_driver.log'), mode='w'),
        log_cfg.default_level
    )
    ntf.add_subscribers(
        (LogNotification, d, cfg.query('notification', 'log_notification')),
        (EmailNotification, d, cfg.query('notification', 'email_notification'))
    )

    log_cfg.switch_formatter(log_cfg.blank_formatter)


def _process_dataset(d):
    """
    :param d: Name of a dataset (not a full path!) to process
    :return: exit status (9 if stacktrace)
    """
    app_logger = get_logger('client')

    dataset_job_dir = os.path.join(cfg['jobs_dir'], d.name)
    if not os.path.isdir(dataset_job_dir):
        os.makedirs(dataset_job_dir)

    setup_logging(d)

    app_logger.info('\nEdinburgh Genomics Analysis Driver')
    with open(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'version.txt'), 'r') as f:
        app_logger.info('Version ' + f.read() + '\n')
    log_cfg.switch_formatter(log_cfg.default_formatter)

    app_logger.info('Using config file at ' + cfg.config_file)
    app_logger.info('Triggering for dataset: ' + d.name)

    exit_status = 9
    stacktrace = None
    try:
        # TODO: launch a pipeline directly which includes the process trigger step instead of launching the process trigger which launch the pipeline
        # Only process the first new dataset found. Run through Cron, this will result
        # in one new pipeline being kicked off per minute.
        from analysis_driver import driver
        ntf.start_pipeline()
        exit_status = driver.pipeline(d)
        app_logger.info('Done')

    except Exception:
        d.fail()
        import traceback
        log_cfg.switch_formatter(log_cfg.blank_formatter)  # blank formatting for stacktrace
        stacktrace = traceback.format_exc()
        log_cfg.switch_formatter(log_cfg.default_formatter)

    finally:
        ntf.end_pipeline(exit_status, stacktrace)
        return exit_status


def _parse_args():
    p = argparse.ArgumentParser()
    p.add_argument('--debug', action='store_true', help='override pipeline log level to debug')
    group = p.add_mutually_exclusive_group(required=True)
    group.add_argument('--run', action='store_true')
    group.add_argument('--sample', action='store_true')
    p.add_argument('--report', action='store_true', help='report on status of datasets')
    p.add_argument('--report-all', action='store_true', help='report all datasets, including finished ones')
    p.add_argument('--skip', nargs='+', default=[], help='mark a dataset as completed')
    p.add_argument('--reset', nargs='+', default=[], help='unmark a dataset as unprocessed for rerunning')
    p.add_argument('--abort', nargs='+', default=[], help='mark a dataset as aborted')
    p.add_argument(
        '--force',
        nargs='+',
        default=[],
        help='mark a sample for processing, even if below the data threshold'
    )

    return p.parse_args()
