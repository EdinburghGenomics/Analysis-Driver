import os
from collections import defaultdict
from egcg_core.rest_communication import get_documents
from egcg_core.app_logging import AppLogger
from egcg_core.util import query_dict
from egcg_core.config import cfg
from analysis_driver.dataset import RunDataset, SampleDataset, ProjectDataset
from egcg_core.constants import DATASET_NEW, DATASET_READY, DATASET_FORCE_READY, DATASET_PROCESSING,\
    DATASET_PROCESSED_SUCCESS, DATASET_PROCESSED_FAIL, DATASET_ABORTED, DATASET_REPROCESS, DATASET_DELETED,\
    DATASET_RESUME


class DatasetScanner(AppLogger):
    type = None
    endpoint = None
    item_id = None

    status_hidden = (DATASET_PROCESSED_SUCCESS, DATASET_PROCESSED_FAIL, DATASET_ABORTED, DATASET_DELETED)
    status_visible = (DATASET_NEW, DATASET_REPROCESS, DATASET_RESUME,
                      DATASET_READY, DATASET_FORCE_READY, DATASET_PROCESSING)

    def __init__(self):
        self.input_dir = cfg.query('%s.input_dir' % self.type)
        self.__triggerignore = None

    def get_dataset(self, *args, **kwargs):
        raise NotImplementedError

    def report(self, all_datasets=False):
        out = [
            '========= %s report =========' % self.__class__.__name__,
            'dataset location: ' + self.input_dir
        ]
        statuses = self.status_visible
        if all_datasets:
            statuses += self.status_hidden
        scan = self.scan_datasets(*statuses)

        for status in statuses:
            datasets = [d.report() for d in scan.get(status, [])]
            if datasets:
                out.append('=== ' + status + ' ===')
                out.append('\n'.join(datasets))

        out.append('_' * 42)
        print('\n'.join(out))

    def _get_dataset_records_for_statuses(self, statuses):
        self.debug('Querying Rest API for status %s', ', '.join(statuses))
        if DATASET_NEW in statuses:
            statuses = list(statuses)
            statuses.remove(DATASET_NEW)
            statuses.append(None)
        if len(statuses) > 1:
            where = {'$or': [{'aggregated.most_recent_proc.status': status} for status in statuses]}
        else:
            where = {'aggregated.most_recent_proc.status': statuses[0]}
        return [
            d for d in get_documents(self.endpoint, where=where, all_pages=True, quiet=True, max_results=100)
            if d[self.item_id] not in self._triggerignore
        ]

    def _get_datasets_for_statuses(self, statuses):
        self.debug('Creating Datasets for status %s', ', '.join(statuses))
        return [
            self.get_dataset(d[self.item_id], query_dict(d, 'aggregated.most_recent_proc'))
            for d in self._get_dataset_records_for_statuses(statuses)
        ]

    def scan_datasets(self, *rest_api_statuses):
        datasets = defaultdict(list)
        self.debug('Scanning for datasets with statuses %s', ', '.join(rest_api_statuses))
        for d in self._get_datasets_for_statuses(rest_api_statuses):
            if d and d.name not in [d.name for d in datasets[d.dataset_status]]:
                datasets[d.dataset_status].append(d)
        for k in datasets:
            datasets[k].sort()

        return datasets

    @property
    def _triggerignore(self):
        if self.__triggerignore is None:
            triggerignore = os.path.join(self.input_dir, '.triggerignore')

            ignorables = []
            if os.path.isfile(triggerignore):
                with open(triggerignore, 'r') as f:
                    for p in f.readlines():
                        if not p.startswith('#'):
                            ignorables.append(p.rstrip('\n'))
            self.debug('Ignoring %s datasets', len(ignorables))
            self.__triggerignore = ignorables
        return self.__triggerignore


class RunScanner(DatasetScanner):
    type = 'run'
    endpoint = 'runs'
    item_id = 'run_id'
    expected_bcl_subdirs = ('RunInfo.xml', 'Data')

    def get_dataset(self, name, most_recent_proc=None):
        dataset_path = os.path.join(self.input_dir, name)
        if os.path.exists(dataset_path):
            return RunDataset(name, most_recent_proc=most_recent_proc)

    def _get_dataset_records_for_statuses(self, statuses):
        rest_api_datasets = super()._get_dataset_records_for_statuses(statuses)
        dataset_names = [d.get(self.item_id) for d in rest_api_datasets]
        if DATASET_NEW in statuses or DATASET_READY in statuses:
            self.debug('Scanning disk for datasets with status %s', ','.join(statuses))
            for d in self._datasets_on_disk():
                if d not in dataset_names:
                    rest_api_datasets.append({self.item_id: d})
        return rest_api_datasets

    def _datasets_on_disk(self):
        return [
            d for d in os.listdir(self.input_dir)
            if self._is_valid_dataset(d) and d not in self._triggerignore
        ]

    def _is_valid_dataset(self, dataset):
        d = os.path.join(self.input_dir, dataset)
        if os.path.isdir(d) and not d.startswith('.'):
            observed = os.listdir(d)
            return all([subdir in observed for subdir in self.expected_bcl_subdirs])
        else:
            return False


class SampleScanner(DatasetScanner):
    type = 'sample'
    endpoint = 'samples'
    item_id = 'sample_id'

    def get_dataset(self, name, most_recent_proc=None):
        return SampleDataset(name, most_recent_proc)


class ProjectScanner(DatasetScanner):
    type = 'project'
    endpoint = 'projects'
    item_id = 'project_id'

    def get_dataset(self, name, most_recent_proc=None):
        return ProjectDataset(name, most_recent_proc)
