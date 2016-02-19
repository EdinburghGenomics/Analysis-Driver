__author__ = 'mwham'
import sys
import os
from os.path import join
import argparse
import logging
import requests

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from analysis_driver.exceptions import AnalysisDriverError
from analysis_driver.config import default as cfg, logging_default as log_cfg
from analysis_driver.app_logging import get_logger
from analysis_driver import rest_communication
from analysis_driver import executor

log_cfg.default_level = logging.DEBUG
log_cfg.add_handler('stdout', logging.StreamHandler(stream=sys.stdout), logging.DEBUG)

app_logger = get_logger('delete_data')


def delete_raw(args):
    non_deleted_runs = query_api(
        'runs',
        'embedded={"run_elements":1,"analysis_driver_procs":1}',
        'aggregate=True'
    )

    deletable_runs = []
    for r in non_deleted_runs:
        if 'not reviewed' not in r.get('review_statuses', ['not reviewed']):
            if r.get('analysis_driver_procs', [{}])[-1].get('status') in ('finished', 'aborted'):
                deletable_runs.append(r)

    app_logger.debug('Found %s runs for deletion' % len(deletable_runs))

    for run in deletable_runs:
        run_id = run['run_id']
        raw_data = join(cfg['input_dir'], run_id)
        deletable_data = join(os.path.dirname(raw_data), '.' + run_id + '.deleting')
        archived_data = join(cfg['achive_dir'], run_id)
        app_logger.debug('Creating deletion dir: ' + deletable_data)

        _execute('mkdir -p ' + deletable_data, args.dry_run)

        deletable_sub_dirs = ('Data', 'Logs', 'Thumbnail_Images')
        for d in deletable_sub_dirs:
            from_d = join(raw_data, d)
            to_d = join(deletable_data, d)

            _execute('mv %s %s' % (from_d, to_d), args.dry_run)

        if not args.dry_run:
            observed = sorted(os.listdir(deletable_data))
            expected = sorted(deletable_sub_dirs)
            assert observed == expected, 'Unexpected deletable sub dirs: ' + str(observed)

        app_logger.debug('Removing deletion dir')
        _execute('rm -rf ' + deletable_data, args.dry_run)

        app_logger.debug('Updating dataset status')
        rest_communication.patch_entry(
            'analysis_driver_procs',
            {'status': 'deleted'},
            proc_id=run['analysis_driver_procs'][-1]['proc_id']
        )

        _execute('mv %s %s' % (raw_data, archived_data), args.dry_run)
        app_logger.debug('Done')


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--dry_run', action='store_true')
    args = p.parse_args()

    cfg.merge(cfg['run'])
    delete_raw(args)


def _execute(cmd, dry_run):
    if dry_run:
        app_logger.info('Running: ' + cmd)
    else:
        status = executor.execute([cmd]).join()
        if status:
            raise AnalysisDriverError('Command failed: ' + cmd)


def query_api(endpoint, *queries):  # TODO: this should go through rest_communication
    query = cfg['rest_api']['url'].rstrip('/') + '/' + endpoint
    if queries:
        query += '?' + '&'.join(queries)
    return requests.get(query).json()['data']


if __name__ == '__main__':
    main()
