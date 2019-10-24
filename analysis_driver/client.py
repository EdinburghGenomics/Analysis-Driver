import os
import sys
import logging
import argparse
import signal
import traceback
from egcg_core import rest_communication
from egcg_core.executor import stop_running_jobs
from egcg_core.app_logging import logging_default as log_cfg
from analysis_driver import exceptions
from analysis_driver.config import default as cfg, load_config
from analysis_driver.dataset_scanner import RunScanner, SampleScanner, ProjectScanner, DATASET_READY,\
    DATASET_FORCE_READY, DATASET_NEW, DATASET_REPROCESS, DATASET_RESUME

app_logger = log_cfg.get_logger('client')


def main(argv=None):
    args = _parse_args(argv)

    load_config()

    log_cfg.set_log_level(logging.DEBUG)
    log_cfg.cfg = cfg.get('logging', {})
    log_cfg.configure_handlers_from_config()

    if args.run:
        scanner = RunScanner()
    elif args.sample:
        scanner = SampleScanner()
    elif args.project:
        scanner = ProjectScanner()
    else:
        raise exceptions.AnalysisDriverError('Invalid arguments: %s' % args)

    if any([args.abort, args.skip, args.reset, args.resume, args.force, args.report, args.report_all, args.stop]):
        log_cfg.add_stdout_handler()

        for d in args.abort:
            scanner.get_dataset(d).abort()
        for d in args.skip:
            scanner.get_dataset(d).skip()
        for d in args.reset:
            scanner.get_dataset(d).reset()
        for d in args.resume:
            scanner.get_dataset(d).resume()
        for d in args.force:
            scanner.get_dataset(d).force()
        for d in args.stop:
            scanner.get_dataset(d).terminate()
        for d in args.soft_stop:
            scanner.get_dataset(d).soft_terminate()

        if args.report:
            scanner.report()
        elif args.report_all:
            scanner.report(all_datasets=True)
        return 0

    processable_statuses = (DATASET_FORCE_READY, DATASET_RESUME, DATASET_READY)
    datasets = scanner.scan_datasets(DATASET_NEW, DATASET_REPROCESS, *processable_statuses)
    ready_datasets = []
    for s in processable_statuses:
        ready_datasets += datasets.get(s, [])

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
        log_cfg.add_handler(logging.FileHandler(repo_log, mode='a'), level=logging.INFO)

    log_cfg.add_handler(
        logging.FileHandler(os.path.join(cfg['jobs_dir'], d.name, 'analysis_driver.log'), mode='w'),
        level=logging.INFO
    )
    log_cfg.add_handler(
        logging.FileHandler(os.path.join(cfg['jobs_dir'], d.name, 'analysis_driver_debug.log'), mode='w'),
        level=logging.DEBUG
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
    with open(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'version.txt'), 'r') as f:
        app_logger.info('\nEdinburgh Genomics Analysis Driver %s\n', f.read())

    log_cfg.set_formatter(log_cfg.default_formatter)

    app_logger.info('Using config file at ' + cfg.config_file)
    app_logger.info('Triggering pipeline %s for dataset %s', d.pipeline.name, d.name)
    app_logger.info('Job dir: %s', dataset_job_dir)

    def _handle_exception(exception):
        app_logger.critical('Encountered a %s exception: %s', exception.__class__.__name__, str(exception))
        etype, value, tb = sys.exc_info()
        if tb:
            stacktrace = ''.join(traceback.format_exception(etype, value, tb))
            app_logger.info('Stacktrace below:\n' + stacktrace)
            d.ntf.crash_report(stacktrace)
        _handle_termination(9)

    def _sigterm_handler(sig, frame):
        app_logger.info('Received signal %s in call stack:\n%s', sig, ''.join(traceback.format_stack(frame)))
        _handle_termination(sig)

    def _handle_termination(sig):
        stop_running_jobs()
        d.fail(sig)
        sys.exit(sig)

    # SIGUSR1 is used by Luigi to stop submitting new jobs.
    # make sure here that SIGUSR1 is caught even if luigi is not running.
    signal.signal(signal.SIGUSR1, _sigterm_handler)
    # SIGUSR2 is used by Analysis Driver to know that another driver process has asked it to terminate.
    signal.signal(signal.SIGUSR2, _sigterm_handler)
    # SIGTERM is used by Analysis Driver to know that a manual process has asked it to be terminated.
    signal.signal(signal.SIGTERM, _sigterm_handler)
    exit_status = 9
    try:
        from analysis_driver import pipelines
        exit_status = pipelines.pipeline(d)
        app_logger.info('Done')

        if exit_status == 0:
            d.succeed()
        else:
            d.fail(exit_status)
        app_logger.info('Finished with exit status ' + str(exit_status))

    except exceptions.SequencingRunError as e:
        app_logger.info('Bad sequencing run: %s. Aborting this dataset', e)
        # Aborting for aborted run is the right thing to do so exit status should be 0
        exit_status = 0
        # Initiate the review as this was not done during the pipeline
        rest_communication.post_entry(
            'actions',
            {'action_type': 'automatic_run_review', 'run_id': d.name},
            use_data=True
        )
        d.abort()
        d.ntf.end_pipeline(exit_status)

    except Exception as e:
        _handle_exception(e)

    finally:
        return exit_status


def _parse_args(argv=None):
    p = argparse.ArgumentParser()
    group = p.add_mutually_exclusive_group(required=True)
    group.add_argument('--run', action='store_true')
    group.add_argument('--sample', action='store_true')
    group.add_argument('--project', action='store_true')
    p.add_argument('--report', action='store_true', help='report on status of datasets')
    p.add_argument('--report-all', action='store_true', help='report all datasets, including finished ones')
    p.add_argument('--skip', nargs='+', default=[], help='mark a dataset as completed')
    p.add_argument('--resume', nargs='+', default=[], help='mark a dataset for rerunning from the last completed stage')
    p.add_argument('--reset', nargs='+', default=[], help='mark a dataset for rerunning from the start')
    p.add_argument('--abort', nargs='+', default=[], help='mark a dataset as aborted')
    p.add_argument('--force', nargs='+', default=[], help='mark a sample for processing, even if below the data threshold')
    p.add_argument('--stop', nargs='+', default=[], help='stop a currently processing run/sample')
    p.add_argument('--soft-stop', nargs='+', default=[], help='prevent any new luigi task from starting for a run/sample')

    return p.parse_args(argv)
