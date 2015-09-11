__author__ = 'mwham'
import argparse
import logging
import os.path
import sys
from analysis_driver import process_trigger as proctrigger
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

    if args.report:
        proctrigger.report()
    elif args.report_all:
        proctrigger.report(all_datasets=True)
    elif args.skip:
        proctrigger.skip(args.skip)
    elif args.reset:
        proctrigger.reset(args.reset)
    else:
        for d, s in proctrigger.scan_datasets():
            if s == 'ready':
                proctrigger.setup_run(d)
                _setup_logging(d, args)
                app_logger = app_logging.get_logger('trigger')
                app_logger.info('Using config file at ' + cfg.config_file)
                app_logger.info('Triggering for dataset: ' + d)

                try:
                    proctrigger.trigger(d)
                    # Only process the first new dataset found. Run through Cron, this will result
                    # in one new pipeline being kicked off per minute.
                    proctrigger.trigger(d)
                    app_logger.info('Done')
                    return 0
                except Exception:
                    _log_stacktrace()
                    return 1

        _setup_logging(None, args)
        app_logger.debug('No new datasets found')

    return 0


def _setup_logging(dataset, args):
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

    # dataset-specific FileHandler
    if dataset:
        dataset_log = os.path.join(cfg['jobs_dir'], dataset, 'analysis_driver.log')
        log_cfg.add_handler('file', logging.FileHandler(filename=dataset_log)) 


def _log_stacktrace():
    import traceback
    lines = traceback.format_exc().splitlines()

    # switch all active logging handlers to a blank Formatter
    log_cfg.switch_formatter(logging.Formatter())
    app_logger = app_logging.get_logger('trigger')
    for line in lines:
        app_logger.error(line)

