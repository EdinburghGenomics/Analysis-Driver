import sys
import os
from os.path import join as p_join
from datetime import datetime
import argparse
import logging

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from analysis_driver.exceptions import AnalysisDriverError
from analysis_driver.config import default as cfg, logging_default as log_cfg
from analysis_driver.app_logging import AppLogger
from analysis_driver import rest_communication, executor, clarity, util


class Deleter(AppLogger):
    def __init__(self, work_dir, dry_run=False, deletion_limit=None):
        self.work_dir = work_dir
        self.dry_run = dry_run
        self.deletion_limit = deletion_limit

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

    def deletable_runs(self):
        runs = rest_communication.get_documents(
            'runs',
            depaginate=True,
            max_results=100,
            embedded={'run_elements': 1, 'analysis_driver_procs': 1},
            sort='run_id',
            aggregate=True
        )

        deletable_runs = []
        for r in runs:
            review_statuses = r.get('review_statuses')
            most_recent_proc = r.get('analysis_driver_procs', [{}])[-1]
            if type(review_statuses) is list:
                review_statuses = [s for s in review_statuses if s]
            if review_statuses and 'not reviewed' not in review_statuses:
                if most_recent_proc.get('status') in ('finished', 'aborted'):  # i.e. not 'deleted'
                    deletable_runs.append(r)

        return deletable_runs[:self.deletion_limit]

    def _setup_run_for_deletion(self, run_id, deletion_dir):
        raw_data = p_join(cfg['input_dir'], run_id)
        deletable_data = p_join(deletion_dir, run_id)
        self.debug('Creating deletion dir: ' + deletable_data)
        self._execute('mkdir -p ' + deletable_data)

        for d in self.deletable_sub_dirs:
            from_d = p_join(raw_data, d)
            to_d = p_join(deletable_data, d)
            self._execute('mv %s %s' % (from_d, to_d))

        return os.listdir(deletable_data)

    def setup_runs_for_deletion(self, runs):
        deletion_dir = p_join(
            cfg['input_dir'],
            '.data_deletion_' + datetime.utcnow().strftime('%d_%m_%Y_%H:%M:%S')
        )
        run_ids = [r['run_id'] for r in runs]
        for run in run_ids:
            deletable_dirs = self._setup_run_for_deletion(run, deletion_dir)
            self._compare_lists(
                deletable_dirs,
                self.deletable_sub_dirs,
                'Unexpected deletable sub dirs:'
            )

        self._compare_lists(os.listdir(deletion_dir), run_ids)
        return deletion_dir

    def delete_runs(self, deletion_dir):
        runs_to_delete = os.listdir(deletion_dir)
        self.debug('Removing deletion dir with %s runs' % len(runs_to_delete))
        self._execute('rm -rfv ' + deletion_dir, cluster_execution=True)

    def mark_run_as_deleted(self, run):
        self.debug('Updating dataset status for ' + run['run_id'])
        if not self.dry_run:
            rest_communication.patch_entry(
                'analysis_driver_procs',
                {'status': 'deleted'},
                'proc_id',
                run['analysis_driver_procs'][-1]['proc_id']
            )

    def archive_run(self, run_id):
        run_to_be_archived = p_join(cfg['input_dir'], run_id)
        self.debug('Archiving ' + run_id)
        assert not any([d in self.deletable_sub_dirs for d in os.listdir(run_to_be_archived)])
        self._execute(
            'mv %s %s' % (
                p_join(cfg['input_dir'], run_id),
                p_join(cfg['data_deletion']['raw_archives'], run_id)
            )
        )
        self.debug('Archiving done')

    def delete_data(self):
        deletable_runs = self.deletable_runs()
        self.debug(
            'Found %s runs for deletion: %s' % (len(deletable_runs), [r['run_id'] for r in deletable_runs])
        )
        if self.dry_run or not deletable_runs:
            return 0

        deletion_dir = self.setup_runs_for_deletion(deletable_runs)
        runs_to_delete = os.listdir(deletion_dir)
        self._compare_lists(
            runs_to_delete,
            [run['run_id'] for run in deletable_runs]
        )
        assert all([os.listdir(p_join(cfg['input_dir'], r)) for r in runs_to_delete])

        for run in deletable_runs:
            assert run['run_id'] in runs_to_delete
            self.mark_run_as_deleted(run)
            self.archive_run(run['run_id'])
            assert os.listdir(p_join(cfg['data_deletion']['raw_archives'], run['run_id']))

        self.delete_runs(deletion_dir)


class FastqDeleter(Deleter):
    def __init__(self, work_dir, dry_run=False, deletion_limit=None):
        super().__init__(work_dir, dry_run, deletion_limit)
        self._samples_released_in_lims = None
        self._samples_released_in_app = None
        # TODO: also check whether the pipeline has delivered the data to Aspera

    @property
    def samples_released_in_lims(self):
        if self._samples_released_in_lims is None:
            self._samples_released_in_lims = clarity.get_released_samples()
        return set(self._samples_released_in_lims)

    @property
    def samples_released_in_app(self):
        if self._samples_released_in_app is None:
            self._samples_released_in_app = rest_communication.get_documents(
                'samples',
                where={'delivered': 'yes', 'useable': 'yes'},  # TODO: do we want useable only?
                projection={'sample_id': 1},
                depaginate=True
            )
        return set([s['sample_id'] for s in self._samples_released_in_app])

    @staticmethod
    def find_run_elements_for_sample(sample_id):
        return rest_communication.get_documents('run_elements', where={'sample_id': sample_id})

    @staticmethod
    def find_fastqs_for_run_element(run_element):
        return util.find_fastqs(
            p_join(cfg['data_deletion']['fastqs'], run_element['run_id'], 'fastq'),
            run_element['project_id'],
            run_element['sample_id'],
            lane=run_element['lane']
        )

    def deletable_samples(self):
        return self.samples_released_in_lims & self.samples_released_in_app

    def find_all_deletable_fastqs(self):
        fastqs = []
        for sample_id in self.deletable_samples():
            for e in self.find_run_elements_for_sample(sample_id):
                fastqs.extend(self.find_fastqs_for_run_element(e))

        return fastqs

    def setup_fastqs_for_deletion(self):
        pass

    def delete_data(self):
        deletable_samples = self.deletable_samples()
        samples_to_delete = {}

        n_samples = 0
        n_run_elements = 0
        n_fastqs = 0

        for s in deletable_samples:
            n_samples += 1
            samples_to_delete[s] = {}
            for e in self.find_run_elements_for_sample(s):
                n_run_elements += 1
                fastqs = self.find_fastqs_for_run_element(e)
                n_fastqs += len(e)
                samples_to_delete[s][e['run_element_id']] = fastqs

        print(samples_to_delete)

        #self.debug(
        print(
            'Found %s deletable fastqs from %s samples and %s run elements' % (
                n_fastqs, n_samples, n_run_elements
            )
        )

        if self.dry_run or not n_fastqs:
            return 0

        # deletion_dir = self.setup_fastqs_for_deletion(deletable_fastqs)
        # fastqs_to_delete = os.listd
        # runs_to_delete = os.listdir(deletion_dir)
        # self._compare_lists(
        #     runs_to_delete,
        #     [run['run_id'] for run in deletable_runs]
        # )
        # assert all([os.listdir(p_join(cfg['input_dir'], r)) for r in runs_to_delete])
        #
        # for run in deletable_runs:
        #     assert run['run_id'] in runs_to_delete
        #     self.mark_run_as_deleted(run)
        #     self.archive_run(run['run_id'])
        #     assert os.listdir(p_join(cfg['archive_dir'], run['run_id']))
        #
        # self.delete_runs(deletion_dir)


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--dry_run', action='store_true')
    p.add_argument('--debug', action='store_true')
    p.add_argument('--work_dir', type=str, required=True)
    p.add_argument('--deletion_limit', type=int, default=None)
    args = p.parse_args()

    if args.debug:
        log_cfg.default_level = logging.DEBUG
        log_cfg.add_handler('stdout', logging.StreamHandler(stream=sys.stdout), logging.DEBUG)

    cfg.merge(cfg['run'])
    d = RawDataDeleter(args.work_dir, args.dry_run, args.deletion_limit)
    d.run_deletion()


if __name__ == '__main__':
    main()
