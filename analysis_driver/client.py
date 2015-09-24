__author__ = 'mwham'
import argparse
import logging
import os
import sys
from analysis_driver import dataset_scanner as scanner
from analysis_driver import app_logging
from analysis_driver.config import default as cfg, logging_default as log_cfg
from analysis_driver.notification import default as ntf, LogNotification, EmailNotification
from analysis_driver.exceptions import AnalysisDriverError


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '--debug',
        action='store_true',
        help='override pipeline log level to debug'
    )
    parser.add_argument(
        '--report',
        action='store_true',
        help='don\'t execute anything, report on status of datasets'
    )
    parser.add_argument(
        '--report-all',
        action='store_true',
        help='report all datasets, including finished ones'
    )
    parser.add_argument(
        '--skip',
        nargs='+',
        default=[],
        help='mark a dataset as completed'
    )
    parser.add_argument(
        '--reset',
        nargs='+',
        default=[],
        help='unmark a dataset as complete/in progress for rerunning'
    )
    parser.add_argument(
        '--abort',
        nargs='+',
        default=[],
        help='mark a dataset as aborted'
    )
    args = parser.parse_args()
    _setup_logging(args)

    if args.reset or args.skip or args.abort:
        for d in args.abort:
            scanner.switch_status(d, 'aborted')
        for d in args.skip:
            scanner.switch_status(d, 'complete')
        for d in args.reset:
            scanner.reset(d)
        return 0

    if args.report:
        scanner.report()
    elif args.report_all:
        scanner.report(all_datasets=True)
    else:
        all_datasets = scanner.scan_datasets()
        new_datasets = all_datasets['new, rta complete']
        if cfg.get('intermediate_dir'):
            new_datasets.extend(all_datasets.get('new', []))

        for d in new_datasets:
            _setup_run(d)
            log_cfg.add_handler(
                'dataset',
                logging.FileHandler(
                    filename=os.path.join(cfg['jobs_dir'], d, 'analysis_driver.log'),
                    mode='w'
                )
            )
            ntf.add_subscribers(
                (LogNotification, d, cfg.query('notification', 'log_notification')),
                (EmailNotification, d, cfg.query('notification', 'email_notification'))
            )

            app_logger = app_logging.get_logger('client')
            _log_app_info(app_logger)
            app_logger.info('Using config file at ' + cfg.config_file)
            app_logger.info('Triggering for dataset: ' + d)

            exit_status = 9
            stacktrace = None
            try:
                # Only process the first new dataset found. Run through Cron, this will result
                # in one new pipeline being kicked off per minute.
                from analysis_driver import process_trigger as proctrigger
                ntf.start_pipeline()
                exit_status = proctrigger.trigger(d)
                app_logger.info('Done')

            except Exception:
                import traceback
                log_cfg.switch_formatter(log_cfg.blank_formatter)  # blank formatting for stacktrace
                stacktrace = traceback.format_exc()

            finally:
                ntf.end_pipeline(exit_status, stacktrace)
                return exit_status

        app_logger = app_logging.get_logger('client')
        app_logger.debug('No new datasets found')

    return 0


def _setup_logging(args):
    if args.debug:
        log_cfg.log_level = logging.DEBUG
    else:
        log_cfg.log_level = logging.INFO

    # user-defined handlers
    logging_handlers = cfg.query('logging', 'handlers')
    if logging_handlers:
        for name, info in logging_handlers.items():
            if 'filename' in info:
                handler = logging.FileHandler(filename=info['filename'])
            elif 'stream' in info:
                s = None
                if info['stream'] == 'ext://sys.stdout':
                    s = sys.stdout
                elif info['stream'] == 'ext://sys.stderr':
                    s = sys.stderr
                handler = logging.StreamHandler(stream=s)
            else:
                raise AnalysisDriverError('Invalid logging configuration: %s %s' % name, str(info))
            log_cfg.add_handler(name, handler)


def _setup_run(dataset):
    for d in [os.path.join(cfg['jobs_dir'], dataset), os.path.join(cfg['jobs_dir'], dataset, 'fastq')]:
        try:
            os.makedirs(d)
        except FileExistsError:
            pass


def _log_app_info(logger):
    log = logger.info
    log_cfg.switch_formatter(log_cfg.blank_formatter)
    log('\nEdinburgh Genomics Analysis Driver')
    version_file = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'version.txt')
    with open(version_file, 'r') as f:
        log('Version ' + f.read() + '\n')
    log_cfg.switch_formatter(log_cfg.default_formatter)
