__author__ = 'mwham'
import sys
import os
from os.path import join
import argparse
import logging

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from analysis_driver.exceptions import AnalysisDriverError
from analysis_driver.config import default as cfg, logging_default as log_cfg
from analysis_driver.app_logging import get_logger
from analysis_driver import rest_communication
from analysis_driver import executor

log_cfg.default_level = logging.DEBUG
log_cfg.add_handler('stdout', logging.StreamHandler(stream=sys.stdout), logging.DEBUG)

app_logger = get_logger('delete_data')


def delete_raw(run_id, dry_run):
    raw_data = join(cfg['input_dir'], run_id)
    dir_for_deletion = join(os.path.dirname(raw_data), '.' + run_id + '.deleting')
    app_logger.debug('Creating deletion dir: ' + dir_for_deletion)
    os.mkdir(dir_for_deletion)

    deletable_sub_dirs = ('Data', 'Logs', 'Thumbnail_Images')
    for d in deletable_sub_dirs:
        from_d = join(raw_data, d)
        to_d = join(dir_for_deletion, d)

        _execute('mv %s %s' % (from_d, to_d), dry_run)

    observed = sorted(os.listdir(dir_for_deletion))
    expected = sorted(deletable_sub_dirs)
    assert observed == expected, 'Unexpected deletable sub dirs: ' + str(observed)
    app_logger.debug('Removing deletion dir')
    _execute('rm -rf ' + dir_for_deletion, dry_run)

    app_logger.debug('Updating dataset status')
    rest_communication.patch_entry('runs', {'status': 'deleted'}, run_id=run_id)
    app_logger.debug('Done')


def main():
    p = argparse.ArgumentParser()
    p.add_argument('run_id')
    p.add_argument('--dry_run', action='store_true')
    args = p.parse_args()

    deletable_runs = rest_communication.get_documents('runs', status='complete', reviewed='yes')
    if deletable_runs:
        delete_raw(deletable_runs[0]['run_id'], dry_run=args.dry_run)


def _execute(cmd, dry_run):
    if dry_run:
        print(cmd)
    else:
        status = executor.execute([cmd]).join()
        if status:
            raise AnalysisDriverError('Command failed: ' + cmd)


if __name__ == '__main__':
    main()
