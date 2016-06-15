from os import listdir
from os.path import basename, expanduser, join, isdir
from datetime import datetime
import argparse
import logging

import sys
from os.path import dirname, abspath
sys.path.append(dirname(dirname(abspath(__file__))))
from analysis_driver.exceptions import AnalysisDriverError
from analysis_driver.config import default as cfg
from analysis_driver.app_logging import AppLogger, logging_default as log_cfg
from analysis_driver.constants import ELEMENT_PROJECT_ID, ELEMENT_SAMPLE_INTERNAL_ID, ELEMENT_PROCS,\
    ELEMENT_RUN_NAME, ELEMENT_STATUS, ELEMENT_PROC_ID, ELEMENT_USEABLE, ELEMENT_FASTQS_DELETED,\
    ELEMENT_DELIVERED, ELEMENT_LANE, DATASET_DELETED
from analysis_driver import executor, util
from analysis_driver.external_data import rest_communication, clarity


class Deleter(AppLogger):
    data_dir = ''

    def __init__(self, work_dir, dry_run=False, deletion_limit=None):
        self.work_dir = work_dir
        self.dry_run = dry_run
        self.deletion_limit = deletion_limit
        self.deletion_dir = join(self.data_dir, '.data_deletion_' + self._strnow())

    def delete_dir(self, d):
        self.debug('Removing deletion dir containing: %s', listdir(d))
        self._execute('rm -rfv ' + d, cluster_execution=True)

    def _execute(self, cmd, cluster_execution=False):
        if cluster_execution:
            status = executor.execute(cmd, job_name='data_deletion', working_dir=self.work_dir).join()
        else:
            status = executor.execute(cmd, env='local').join()
        if status:
            raise AnalysisDriverError('Command failed: ' + cmd)

    def _compare_lists(self, observed, expected, error_message='List comparison mismatch:'):
        observed = sorted(observed)
        expected = sorted(expected)
        if observed != expected:
            self.error(error_message)
            self.error('observed: ' + str(observed))
            self.error('expected: ' + str(expected))
            raise AssertionError

    @staticmethod
    def _now():
        return datetime.utcnow()

    @classmethod
    def _strnow(cls):
        return cls._now().strftime('%d_%m_%Y_%H:%M:%S')


class RawDataDeleter(Deleter):
    deletable_sub_dirs = ('Data', 'Logs', 'Thumbnail_Images')
    data_dir = cfg['data_deletion']['raw_data']

    def __init__(self, work_dir, dry_run=False, deletion_limit=None, manual_delete=None):
        super().__init__(work_dir, dry_run, deletion_limit)
        self.archive_dir = cfg['data_deletion']['raw_archives']
        self.list_runs = manual_delete

    def deletable_runs(self):
        runs = rest_communication.get_documents('aggregate/all_runs', paginate=False, quiet=True, sort=ELEMENT_RUN_NAME)
        deletable_runs = []
        for r in runs:
            if (self.list_runs and r[ELEMENT_RUN_NAME] in self.list_runs) or self._run_deletable(r):
                deletable_runs.append(r)
        return deletable_runs[:self.deletion_limit]

    def _setup_run_for_deletion(self, run_id):
        raw_data = join(self.data_dir, run_id)
        deletable_data = join(self.deletion_dir, run_id)
        self.debug('Creating deletion dir: ' + deletable_data)
        self._execute('mkdir -p ' + deletable_data)

        for d in self.deletable_sub_dirs:
            from_d = join(raw_data, d)
            if isdir(from_d):
                to_d = join(deletable_data, d)
                self._execute('mv %s %s' % (from_d, to_d))

        return listdir(deletable_data)  # Data, Thumbnail_Images, etc.

    def setup_runs_for_deletion(self, runs):
        run_ids = [r[ELEMENT_RUN_NAME] for r in runs]
        for run in run_ids:
            deletable_dirs = self._setup_run_for_deletion(run)
            if sorted(self.deletable_sub_dirs) != sorted(deletable_dirs):
                self.warning(
                    'Not all deletable dirs were present for run %s: %s' % (
                        run, [d for d in self.deletable_sub_dirs if d not in deletable_dirs]
                    )
                )

        self._compare_lists(listdir(self.deletion_dir), run_ids)

    def mark_run_as_deleted(self, run):
        self.debug('Updating dataset status for ' + run[ELEMENT_RUN_NAME])
        if not self.dry_run:
            rest_communication.patch_entry(
                ELEMENT_PROCS,
                {ELEMENT_STATUS: DATASET_DELETED},
                ELEMENT_PROC_ID,
                run['most_recent_proc'][ELEMENT_PROC_ID]
            )

    def archive_run(self, run_id):
        run_to_be_archived = join(self.data_dir, run_id)
        self.debug('Archiving ' + run_id)
        assert not any([d in self.deletable_sub_dirs for d in listdir(run_to_be_archived)])
        self._execute('mv %s %s' % (join(self.data_dir, run_id), join(self.archive_dir, run_id)))

    def delete_data(self):
        deletable_runs = self.deletable_runs()
        self.debug(
            'Found %s runs for deletion: %s' % (
                len(deletable_runs), [r[ELEMENT_RUN_NAME] for r in deletable_runs]
            )
        )
        if self.dry_run or not deletable_runs:
            return 0

        self.setup_runs_for_deletion(deletable_runs)
        runs_to_delete = listdir(self.deletion_dir)
        self._compare_lists(
            runs_to_delete,
            [run[ELEMENT_RUN_NAME] for run in deletable_runs]
        )
        assert all([listdir(join(self.data_dir, r)) for r in runs_to_delete])

        for run in deletable_runs:
            assert run[ELEMENT_RUN_NAME] in runs_to_delete
            self.mark_run_as_deleted(run)
            self.archive_run(run[ELEMENT_RUN_NAME])
            assert listdir(join(self.archive_dir, run[ELEMENT_RUN_NAME]))

        self.delete_dir(self.deletion_dir)

    @staticmethod
    def _run_deletable(e):
        review_statuses = e.get('review_statuses')
        if type(review_statuses) is list:
            review_statuses = [s for s in review_statuses if s]
        if review_statuses and 'not reviewed' not in review_statuses:  # all run elements have been reviewed
            if e.get('most_recent_proc', {}).get('status') in ('finished', 'aborted'):  # run is not deleted
                return True
        return False


class FastqDeleter(Deleter):
    data_dir = cfg['data_deletion']['fastqs']

    def __init__(self, work_dir, dry_run=False, deletion_limit=None, project_id=None):
        super().__init__(work_dir, dry_run, deletion_limit)
        self._samples_released_in_lims = None
        self._samples_released_in_app = None
        self.project_id = project_id

    @property
    def samples_released_in_lims(self):
        if self._samples_released_in_lims is None:
            self._samples_released_in_lims = clarity.get_released_samples()
        return set(self._samples_released_in_lims)

    @property
    def samples_released_in_app(self):
        if self._samples_released_in_app is None:
            where = {ELEMENT_DELIVERED: 'yes', ELEMENT_USEABLE: 'yes', ELEMENT_FASTQS_DELETED: 'no'}
            if self.project_id:
                where[ELEMENT_PROJECT_ID] = self.project_id
            self._samples_released_in_app = rest_communication.get_documents(
                'samples',
                quiet=True,
                where=where,  # TODO: do we want useable only?
                projection={ELEMENT_SAMPLE_INTERNAL_ID: 1},
                depaginate=True
            )
        return [s[ELEMENT_SAMPLE_INTERNAL_ID] for s in self._samples_released_in_app]

    def find_fastqs_for_run_element(self, run_element):
        return util.find_fastqs(
            join(self.data_dir, run_element[ELEMENT_RUN_NAME], 'fastq'),
            run_element[ELEMENT_PROJECT_ID],
            run_element[ELEMENT_SAMPLE_INTERNAL_ID],
            lane=run_element[ELEMENT_LANE]
        )

    def deletable_samples(self):
        return sorted(self.samples_released_in_lims & set(self.samples_released_in_app))[:self.deletion_limit]

    def setup_deletion_records(self):
        deletable_samples = self.deletable_samples()
        deletion_records = []

        n_samples = 0
        n_run_elements = 0
        n_fastqs = 0

        for s in deletable_samples:
            n_samples += 1
            run_elements = rest_communication.get_documents(
                'run_elements', quiet=True, where={ELEMENT_SAMPLE_INTERNAL_ID: s}
            )
            for e in run_elements:
                assert e[ELEMENT_SAMPLE_INTERNAL_ID] == s
                fastqs = self.find_fastqs_for_run_element(e)
                if fastqs:
                    n_run_elements += 1
                    n_fastqs += len(fastqs)
                    deletion_records.append(_FastqDeletionRecord(e, fastqs))

        self.debug(
            'Found %s deletable fastqs from %s run elements in %s samples: %s' % (
                n_fastqs, n_run_elements, n_samples, deletion_records
            )
        )
        return deletion_records

    def setup_fastqs_for_deletion(self, deletion_records):
        all_fastqs = []
        for r in deletion_records:
            all_fastqs.extend(self._setup_record_for_deletion(r))

        comparisons = (
            (listdir(self.deletion_dir), [r.run_id for r in deletion_records]),
            (util.find_files(self.deletion_dir, '*', 'fastq', '*'), [r.project_id for r in deletion_records]),
            (util.find_files(self.deletion_dir, '*', 'fastq', '*', '*'), [r.sample_id for r in deletion_records]),
            (util.find_all_fastqs(self.deletion_dir), all_fastqs)
        )
        for expected, observed in comparisons:
            self._compare_lists(set([basename(f) for f in observed]),
                                set([basename(f) for f in expected]))

    def _setup_record_for_deletion(self, del_record):
        deletion_sub_dir = join(
            self.deletion_dir,
            del_record.run_id,
            'fastq',
            del_record.project_id,
            del_record.sample_id
        )
        self._execute('mkdir -p ' + deletion_sub_dir)
        self._execute('mv %s %s' % (' '.join(del_record.fastqs), deletion_sub_dir))
        observed = [basename(f) for f in listdir(deletion_sub_dir)]
        assert all([basename(f) in observed for f in del_record.fastqs])
        return del_record.fastqs

    def mark_sample_as_deleted(self, sample_id):
        self.debug('Updating dataset status for ' + sample_id)
        if not self.dry_run:
            rest_communication.patch_entry(
                'samples', {ELEMENT_FASTQS_DELETED: 'yes'}, ELEMENT_SAMPLE_INTERNAL_ID, sample_id
            )

    def delete_data(self):
        deletion_records = self.setup_deletion_records()
        if self.dry_run or not deletion_records:
            return 0

        self.setup_fastqs_for_deletion(deletion_records)
        for s in set([r.sample_id for r in deletion_records]):
            self.mark_sample_as_deleted(s)

        self.delete_dir(self.deletion_dir)


class DeliveredDataDeleter(Deleter):
    data_dir = cfg['data_deletion']['delivered_data']

    def __init__(self, work_dir, dry_run=False, deletion_limit=None, manual_delete=None, sample_ids=None):
        super().__init__(work_dir, dry_run, deletion_limit)
        if manual_delete is None:
            manual_delete = []
        self.list_samples = manual_delete
        self.limit_samples = sample_ids

    def deletable_samples(self):
        manual_samples = []
        if self.list_samples:
            manual_samples = rest_communication.get_documents(
                'aggregate/samples',
                quiet=True,
                match={'$or': [{'sample_id': s} for s in self.list_samples]}
            )
        return sorted(manual_samples + self._auto_deletable_samples(), key=lambda e: e['sample_id'])

    def _auto_deletable_samples(self):
        samples = []
        sample_records = rest_communication.get_documents(
            'aggregate/samples',
            quiet=True,
            match={'$or': [{'proc_status': 'finished'}, {'proc_status': 'aborted'}]}
        )
        for r in sample_records:
            release_date = clarity.get_sample_release_date(r['sample_id'])
            if release_date and self._old_enough_for_deletion(release_date):
                samples.append(r)
        return samples

    def mark_sample_as_deleted(self, sample_id, proc_id):
        if not proc_id:
            self.warning('No pipeline process found for ' + sample_id)
        else:
            rest_communication.patch_entry(
                'analysis_driver_procs',
                {'status': 'deleted'},
                'proc_id', proc_id
            )

    def setup_samples_for_deletion(self, sample_records):
        for s in sample_records:
            self._setup_sample_for_deletion(s['project_id'], s['sample_id'])

    @classmethod
    def _old_enough_for_deletion(cls, date_run, age_threshold=90):
        year, month, day = date_run.split('-')
        age = cls._now() - datetime(int(year), int(month), int(day))
        return age.days > age_threshold

    def _setup_sample_for_deletion(self, project_id, sample_id):
        release_folders = util.find_files(self.data_dir, project_id, '*', sample_id)
        if len(release_folders) != 1:
            self.warning(
                'Found %s deletable directories for sample %s: %s',
                len(release_folders),
                sample_id,
                release_folders
            )
        else:
            deletable_data = join(self.deletion_dir, project_id)
            self._execute('mkdir -p ' + deletable_data)
            self._execute('mv %s %s' % (release_folders[0], deletable_data))

    def _try_delete_empty_dir(self, d):
        if isdir(d) and not listdir(d):
            self._execute('rm -r ' + d)

    def _cleanup_empty_dirs(self, project_id):
        project_folder = join(self.data_dir, project_id)
        for release_dir in listdir(project_folder):
            self._try_delete_empty_dir(join(project_folder, release_dir))
        self._try_delete_empty_dir(project_folder)

    def delete_data(self):
        deletable_samples = self.deletable_samples()
        if self.limit_samples:
            deletable_samples = [s for s in deletable_samples if s['sample_id'] in self.limit_samples]

        sample_ids = [e[ELEMENT_SAMPLE_INTERNAL_ID] for e in deletable_samples]
        self.debug('Found %s samples for deletion: %s' % (len(deletable_samples), sample_ids))

        if self.dry_run or not deletable_samples:
            return 0

        self.setup_samples_for_deletion(deletable_samples)
        samples_to_delete = []
        for p in listdir(self.deletion_dir):
            samples_to_delete.extend(listdir(join(self.deletion_dir, p)))
        self._compare_lists(samples_to_delete, sample_ids)

        for s in deletable_samples:
            self.mark_sample_as_deleted(s['sample_id'], s.get('most_recent_proc', {}).get('proc_id'))

        self.delete_dir(self.deletion_dir)
        for p in set([e[ELEMENT_PROJECT_ID] for e in deletable_samples]):
            self._cleanup_empty_dirs(p)


class _FastqDeletionRecord:
    def __init__(self, run_element, fastqs):
        self.run_element = run_element
        self.fastqs = fastqs
        self.run_id = run_element[ELEMENT_RUN_NAME]
        self.sample_id = run_element[ELEMENT_SAMPLE_INTERNAL_ID]
        self.project_id = run_element[ELEMENT_PROJECT_ID]
        self.lane = run_element[ELEMENT_LANE]

    def __repr__(self):
        return '%s(%s/%s/%s/%s)' % (
            self.__class__.__name__,
            self.run_id,
            self.project_id,
            self.sample_id,
            [basename(f) for f in self.fastqs]
        )


def main():
    deleters = {
        'raw': RawDataDeleter,
        'fastq': FastqDeleter,
        'delivered_data': DeliveredDataDeleter
    }

    p = argparse.ArgumentParser()
    p.add_argument('deleter', type=str, choices=deleters.keys())
    p.add_argument('--debug', action='store_true')
    p.add_argument('--dry_run', action='store_true')
    p.add_argument('--work_dir', default=expanduser('~'))
    p.add_argument('--deletion_limit', type=int, default=None)
    p.add_argument('--project_id', type=str)
    p.add_argument('--manual_delete', type=str, nargs='+')
    p.add_argument('--sample_ids', nargs='+', default=[])
    args = p.parse_args()

    if args.__dict__.pop('debug', False):
        log_level = logging.DEBUG
    else:
        log_level = logging.WARNING

    log_cfg.set_log_level(log_level)
    log_cfg.add_handler(logging.StreamHandler(stream=sys.stdout), log_level)

    deleter_type = args.__dict__.pop('deleter')
    deleter_args = dict([(k, v) for k, v in args.__dict__.items() if v])

    d = deleters[deleter_type](**deleter_args)
    d.delete_data()


if __name__ == '__main__':
    main()
