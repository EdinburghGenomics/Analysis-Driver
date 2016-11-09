import os
import sys
import logging
import argparse
import signal
from egcg_core.executor import stop_running_jobs
from egcg_core.app_logging import logging_default as log_cfg
from analysis_driver import exceptions
from analysis_driver.config import default as cfg, load_config
from analysis_driver.notification import default as ntf, LogNotification, EmailNotification, AsanaNotification
from analysis_driver.dataset_scanner import RunScanner, SampleScanner, DATASET_READY, DATASET_FORCE_READY, DATASET_NEW, DATASET_REPROCESS

app_logger = log_cfg.get_logger('client')


def main():
    args = _parse_args()

    load_config()

    if args.debug:
        log_cfg.log_level = logging.DEBUG

    log_cfg.cfg = cfg.get('logging', {})
    log_cfg.configure_handlers_from_config()

    if args.run:
        if 'run' in cfg:
            cfg.merge(cfg['run'])
        scanner = RunScanner(cfg)
    else:
        assert args.sample
        if 'sample' in cfg:
            cfg.merge(cfg['sample'])
        scanner = SampleScanner(cfg)

    if any([args.abort, args.skip, args.reset, args.force, args.report, args.report_all, args.stop]):
        for d in args.abort:
            scanner.get_dataset(d).abort()
        for d in args.skip:
            dataset = scanner.get_dataset(d)
            dataset.reset()
            dataset.start()
            dataset.succeed(quiet=True)
        for d in args.reset:
            scanner.get_dataset(d).reset()
        for d in args.force:
            scanner.get_dataset(d).force()
        for d in args.stop:
            pid = scanner.get_dataset(d).most_recent_proc.get('pid')
            if pid:
                os.kill(pid, 10)

        if args.report:
            scanner.report()
        elif args.report_all:
            scanner.report(all_datasets=True)
        return 0

    datasets = scanner.scan_datasets(DATASET_NEW, DATASET_REPROCESS, DATASET_READY, DATASET_FORCE_READY)
    ready_datasets = datasets.get(DATASET_FORCE_READY, []) + datasets.get(DATASET_READY, [])
    if not ready_datasets:
        return 0
    else:
        # Only process the first new dataset found. Run through Cron, this will result in one new pipeline
        # being kicked off per minute.
        return _process_dataset(ready_datasets[0])


def setup_dataset_logging(d):
    log_repo = cfg.query('logging', 'repo')
    if log_repo:
        repo_log = os.path.join(log_repo, d.name + '.log')
        log_cfg.add_handler(logging.FileHandler(filename=repo_log, mode='a'))

    job_dir_log = os.path.join(cfg['jobs_dir'], d.name, 'analysis_driver.log')
    log_cfg.add_handler(
        logging.FileHandler(filename=job_dir_log, mode='w')
    )


def _process_dataset(d):
    """
    :param Dataset d: Run or Sample to process
    :return: exit status (9 if stacktrace)
    """
    dataset_job_dir = os.path.join(cfg['jobs_dir'], d.name)
    if not os.path.isdir(dataset_job_dir):
        os.makedirs(dataset_job_dir)

    setup_dataset_logging(d)

    log_cfg.set_formatter(log_cfg.blank_formatter)
    app_logger.info('\nEdinburgh Genomics Analysis Driver')
    with open(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'version.txt'), 'r') as f:
        app_logger.info('Version ' + f.read() + '\n')
    log_cfg.set_formatter(log_cfg.default_formatter)

    app_logger.info('Using config file at ' + cfg.config_file)
    app_logger.info('Triggering for dataset: ' + d.name)

    ntf.add_subscribers(
        (LogNotification, d, cfg.query('notification', 'log_notification')),
        (EmailNotification, d, cfg.query('notification', 'email_notification')),
        (AsanaNotification, d, cfg.query('notification', 'asana'))
    )

    def _handle_exception(exception):
        app_logger.critical('Encountered a %s exception: %s', exception.__class__.__name__, str(exception))
        _handle_termination(9)

    def _sigterm_handler(sig, frame):
        app_logger.info('Received signal %s. Cleaning up running jobs', sig)
        _handle_termination(sig)

    def _handle_termination(_exit_status):
        stop_running_jobs()
        d.fail(_exit_status)

        if any(sys.exc_info()):
            import traceback
            stacktrace = traceback.format_exc()
            app_logger.info('Stack trace below:\n' + stacktrace)
            ntf.crash_report(stacktrace)
        else:
            app_logger.info('No stacktrace to log')

        sys.exit(_exit_status)

    signal.signal(10, _sigterm_handler)
    signal.signal(15, _sigterm_handler)
    exit_status = 9
    try:
        from analysis_driver import driver
        d.start()
        exit_status = driver.pipeline(d)
        app_logger.info('Done')

    except exceptions.SequencingRunError as e:
        app_logger.info('Bad sequencing run: %s. Aborting this dataset' % str(e))
        exit_status = 2  # TODO: we should send a notification of the run status found
        d.abort()

    except Exception as e:
        _handle_exception(e)

    else:
        if exit_status == 0:
            d.succeed()
        else:
            d.fail(exit_status)
        app_logger.info('Finished with exit status ' + str(exit_status))

    finally:
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
    p.add_argument('--stop', nargs='+', default=[], help='stop a currently processing run/sample')

    return p.parse_args()
