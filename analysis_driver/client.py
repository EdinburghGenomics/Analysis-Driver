__author__ = 'mwham'
import argparse
import logging
import os.path
import sys
from analysis_driver import dataset_scanner as scanner
from analysis_driver import app_logging
from analysis_driver.config import default as cfg, logging_default as log_cfg
from analysis_driver.exceptions import AnalysisDriverError


def main():
    parser = argparse.ArgumentParser()
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

    _setup_logging(args)

    if args.report:
        scanner.report()
    elif args.report_all:
        scanner.report(all_datasets=True)
    elif args.skip:
        scanner.skip(args.skip)
    elif args.reset:
        scanner.reset(args.reset)
    else:
        if cfg.get('intermediate_dir'):
            use_int_dir = True
            required_status = 'new'
        else:
            use_int_dir = False
            required_status = 'new, rta complete'
        for d, s in scanner.scan_datasets():
            if required_status in s:
                app_logger = app_logging.get_logger('client')
                log_cfg.add_handler(
                    'dataset',
                    logging.FileHandler(filename=os.path.join(cfg['jobs_dir'], d, 'analysis_driver.log'))
                )

                try:
                    from analysis_driver import process_trigger as proctrigger
                    proctrigger.setup_run(d)
                    app_logger.info('Using config file at ' + cfg.config_file)
                    app_logger.info('Triggering for dataset: ' + d)
                    # Only process the first new dataset found. Run through Cron, this will result
                    # in one new pipeline being kicked off per minute.
                    proctrigger.trigger(d, use_int_dir)
                    app_logger.info('Done')
                    return 0

                except Exception:
                    _log_stacktrace(app_logger)
                    return 1

        app_logger = app_logging.get_logger('client')
        app_logger.debug('No new datasets found')

    return 0


def _setup_logging(args):
    if args.debug:
        log_cfg.log_level = logging.DEBUG
    else:
        log_cfg.log_level = logging.INFO

    # user-defined handlers
    if cfg['logging']['handlers']:
        for name, info in cfg['logging']['handlers'].items():
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


def _log_stacktrace(logger):
    import traceback
    lines = traceback.format_exc().splitlines()

    # switch all active logging handlers to a blank Formatter
    log_cfg.switch_formatter(logging.Formatter())
    for line in lines:
        logger.error(line)
