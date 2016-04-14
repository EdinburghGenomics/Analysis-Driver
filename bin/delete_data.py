from os import listdir
from os.path import basename, expanduser, join as p_join, isdir
from datetime import datetime
import argparse
import logging

import sys
from os.path import dirname, abspath
sys.path.append(dirname(dirname(abspath(__file__))))
from analysis_driver.exceptions import AnalysisDriverError
from analysis_driver.config import default as cfg
from analysis_driver.app_logging import AppLogger, logging_default as log_cfg
from analysis_driver.constants import ELEMENT_PROJECT_ID, ELEMENT_SAMPLE_INTERNAL_ID, ELEMENT_RUN_ELEMENTS,\
    ELEMENT_PROCS, ELEMENT_RUN_NAME, ELEMENT_STATUS, ELEMENT_PROC_ID, ELEMENT_USEABLE, ELEMENT_FASTQS_DELETED,\
    ELEMENT_DELIVERED, ELEMENT_LANE, DATASET_DELETED
from analysis_driver import rest_communication, executor, clarity, util


class Deleter(AppLogger):
    def __init__(self, work_dir, dry_run=False, deletion_limit=None):
        self.work_dir = work_dir
        self.dry_run = dry_run
        self.deletion_limit = deletion_limit

    def delete_dir(self, d):
        self.debug('Removing deletion dir containing: %s', listdir(d))
        self._execute('rm -rfv ' + d, cluster_execution=True)

    def _execute(self, cmd, cluster_execution=False):
        if cluster_execution:
            status = executor.execute([cmd], job_name='data_deletion', working_dir=self.work_dir).join()
        else:
            status = executor.execute([cmd], 'local').join()
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


class RawDataDeleter(Deleter):
    deletable_sub_dirs = ('Data', 'Logs', 'Thumbnail_Images')

    def __init__(self, work_dir, dry_run=False, deletion_limit=None):
        super().__init__(work_dir, dry_run, deletion_limit)
        self.data_dir = cfg['data_deletion']['raw_data']
        self.archive_dir = cfg['data_deletion']['raw_archives']

    def deletable_runs(self):
        runs = rest_communication.get_documents(
            'runs',
            depaginate=True,
            max_results=100,
            embedded={ELEMENT_RUN_ELEMENTS: 1, ELEMENT_PROCS: 1},
            sort=ELEMENT_RUN_NAME,
            aggregate=True
        )

        deletable_runs = []
        for r in runs:
            review_statuses = r.get('review_statuses')
            most_recent_proc = r.get(ELEMENT_PROCS, [{}])[-1]
            if type(review_statuses) is list:
                review_statuses = [s for s in review_statuses if s]
            if review_statuses and 'not reviewed' not in review_statuses:
                if most_recent_proc.get(ELEMENT_STATUS) in ('finished', 'aborted'):  # i.e. not 'deleted'
                    deletable_runs.append(r)

        return deletable_runs[:self.deletion_limit]

    def _setup_run_for_deletion(self, run_id, deletion_dir):
        raw_data = p_join(self.data_dir, run_id)
        deletable_data = p_join(deletion_dir, run_id)
        self.debug('Creating deletion dir: ' + deletable_data)
        self._execute('mkdir -p ' + deletable_data)

        for d in self.deletable_sub_dirs:
            from_d = p_join(raw_data, d)
            if isdir(from_d):
                to_d = p_join(deletable_data, d)
                self._execute('mv %s %s' % (from_d, to_d))

        return listdir(deletable_data)  # Data, Thumbnail_Images, etc.

    def setup_runs_for_deletion(self, runs):
        deletion_dir = p_join(
            self.data_dir,
            '.data_deletion_' + datetime.utcnow().strftime('%d_%m_%Y_%H:%M:%S')
        )
        run_ids = [r[ELEMENT_RUN_NAME] for r in runs]
        for run in run_ids:
            deletable_dirs = self._setup_run_for_deletion(run, deletion_dir)
            if sorted(self.deletable_sub_dirs) != sorted(deletable_dirs):
                self.warning(
                    'Not all deletable dirs were present for run %s: %s' % (
                        run, [d for d in self.deletable_sub_dirs if d not in deletable_dirs]
                    )
                )

        self._compare_lists(listdir(deletion_dir), run_ids)
        return deletion_dir

    def mark_run_as_deleted(self, run):
        self.debug('Updating dataset status for ' + run[ELEMENT_RUN_NAME])
        if not self.dry_run:
            rest_communication.patch_entry(
                ELEMENT_PROCS,
                {ELEMENT_STATUS: DATASET_DELETED},
                ELEMENT_PROC_ID,
                run[ELEMENT_PROCS][-1][ELEMENT_PROC_ID]
            )

    def archive_run(self, run_id):
        run_to_be_archived = p_join(self.data_dir, run_id)
        self.debug('Archiving ' + run_id)
        assert not any([d in self.deletable_sub_dirs for d in listdir(run_to_be_archived)])
        self._execute(
            'mv %s %s' % (
                p_join(self.data_dir, run_id),
                p_join(self.archive_dir, run_id)
            )
        )

    def delete_data(self):
        deletable_runs = self.deletable_runs()
        self.debug(
            'Found %s runs for deletion: %s' % (
                len(deletable_runs), [r[ELEMENT_RUN_NAME] for r in deletable_runs]
            )
        )
        if self.dry_run or not deletable_runs:
            return 0

        deletion_dir = self.setup_runs_for_deletion(deletable_runs)
        runs_to_delete = listdir(deletion_dir)
        self._compare_lists(
            runs_to_delete,
            [run[ELEMENT_RUN_NAME] for run in deletable_runs]
        )
        assert all([listdir(p_join(self.data_dir, r)) for r in runs_to_delete])

        for run in deletable_runs:
            assert run[ELEMENT_RUN_NAME] in runs_to_delete
            self.mark_run_as_deleted(run)
            self.archive_run(run[ELEMENT_RUN_NAME])
            assert listdir(p_join(self.archive_dir, run[ELEMENT_RUN_NAME]))

        self.delete_dir(deletion_dir)


class FastqDeleter(Deleter):
    def __init__(self, work_dir, dry_run=False, deletion_limit=None, project_id=None):
        super().__init__(work_dir, dry_run, deletion_limit)
        self.data_dir = cfg['data_deletion']['fastqs']
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
                where=where,  # TODO: do we want useable only?
                projection={ELEMENT_SAMPLE_INTERNAL_ID: 1},
                depaginate=True
            )
        return [s[ELEMENT_SAMPLE_INTERNAL_ID] for s in self._samples_released_in_app]

    def find_fastqs_for_run_element(self, run_element):
        return util.find_fastqs(
            p_join(self.data_dir, run_element[ELEMENT_RUN_NAME], 'fastq'),
            run_element[ELEMENT_PROJECT_ID],
            run_element[ELEMENT_SAMPLE_INTERNAL_ID],
            lane=run_element[ELEMENT_LANE]
        )

    def setup_deletion_records(self):
        deletable_samples = sorted(self.samples_released_in_lims & set(self.samples_released_in_app))[:self.deletion_limit]
        deletion_records = []

        n_samples = 0
        n_run_elements = 0
        n_fastqs = 0

        for s in deletable_samples:
            n_samples += 1
            for e in rest_communication.get_documents('run_elements', where={ELEMENT_SAMPLE_INTERNAL_ID: s}):
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
        deletion_dir = p_join(
            self.data_dir,
            '.data_deletion_' + datetime.utcnow().strftime('%d_%m_%Y_%H:%M:%S')
        )
        all_fastqs = []
        for r in deletion_records:
            all_fastqs.extend(self._setup_record_for_deletion(deletion_dir, r))

        comparisons = (
            (listdir(deletion_dir), [r.run_id for r in deletion_records]),
            (util.find_files(deletion_dir, '*', 'fastq', '*'), [r.project_id for r in deletion_records]),
            (util.find_files(deletion_dir, '*', 'fastq', '*', '*'), [r.sample_id for r in deletion_records]),
            (util.find_all_fastqs(deletion_dir), all_fastqs)
        )
        for expected, observed in comparisons:
            self._compare_lists(set([basename(f) for f in observed]),
                                set([basename(f) for f in expected]))

        return deletion_dir

    def _setup_record_for_deletion(self, deletion_dir, del_record):
        deletion_sub_dir = p_join(
            deletion_dir,
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

    def mark_sample_fastqs_as_deleted(self, sample_id):
        self.debug('Updating dataset status for ' + sample_id)
        if not self.dry_run:
            rest_communication.patch_entry(
                'samples', {ELEMENT_FASTQS_DELETED: 'yes'}, ELEMENT_SAMPLE_INTERNAL_ID, sample_id
            )

    def delete_data(self):
        deletion_records = self.setup_deletion_records()
        if self.dry_run or not deletion_records:
            return 0

        deletion_dir = self.setup_fastqs_for_deletion(deletion_records)
        for s in set([r.sample_id for r in deletion_records]):
            self.mark_sample_fastqs_as_deleted(s)

        self.delete_dir(deletion_dir)


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
        'fastq': FastqDeleter  # ,
        # 'delivered_data': DeliveredDataDeleter
    }

    p = argparse.ArgumentParser()
    p.add_argument('deleter', type=str, choices=deleters.keys())
    p.add_argument('--debug', action='store_true')
    p.add_argument('--dry_run', action='store_true')
    p.add_argument('--work_dir', default=expanduser('~'))
    p.add_argument('--deletion_limit', type=int, default=None)
    p.add_argument('--project_id', type=str)
    args = p.parse_args()

    if args.__dict__.pop('debug', False):
        log_cfg.default_level = logging.DEBUG
        log_cfg.add_handler(logging.StreamHandler(stream=sys.stdout), logging.DEBUG)

    deleter_type = args.__dict__.pop('deleter')
    deleter_args = dict([(k, v) for k, v in args.__dict__.items() if v])

    d = deleters[deleter_type](**deleter_args)
    d.delete_data()


if __name__ == '__main__':
    main()
