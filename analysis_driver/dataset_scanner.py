import os
from collections import defaultdict
from analysis_driver.app_logging import AppLogger
from analysis_driver.external_data import rest_communication
from analysis_driver.dataset import RunDataset, SampleDataset
from analysis_driver.external_data.clarity import get_list_of_samples, sanitize_user_id
from analysis_driver.constants import DATASET_NEW, DATASET_READY, DATASET_FORCE_READY, DATASET_PROCESSING,\
    DATASET_PROCESSED_SUCCESS, DATASET_PROCESSED_FAIL, DATASET_ABORTED, DATASET_REPROCESS, DATASET_DELETED


class DatasetScanner(AppLogger):
    endpoint = None
    item_id = None

    status_visible = (DATASET_NEW, DATASET_REPROCESS, DATASET_READY, DATASET_FORCE_READY, DATASET_PROCESSING)
    status_hidden = (DATASET_PROCESSED_SUCCESS, DATASET_PROCESSED_FAIL, DATASET_ABORTED, DATASET_DELETED)

    def __init__(self, config):
        self.input_dir = config.get('input_dir')
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
            datasets = [str(d) for d in scan.get(status, [])]
            if datasets:
                out.append('=== ' + status + ' ===')
                out.append('\n'.join(datasets))

        out.append('_' * 42)
        print('\n'.join(out))

    def _get_dataset_records_for_status(self, status):
        self.debug('Querying Rest API for status %s', status)
        if status == DATASET_NEW:
            status = None
        return [
            d for d in rest_communication.get_documents(self.endpoint, match={'proc_status': status}, paginate=False)
            if d[self.item_id] not in self._triggerignore
            ]

    def _get_datasets_for_status(self, status):
        self.debug('Creating Datasets for status %s', status)
        return [
            self.get_dataset(d[self.item_id], d.get('most_recent_proc'))
            for d in self._get_dataset_records_for_status(status)
        ]

    def scan_datasets(self, *rest_api_statuses):
        datasets = defaultdict(list)
        for s in rest_api_statuses:
            self.debug('Scanning for datasets with status %s', s)
            for d in self._get_datasets_for_status(s):
                if d.name not in [d.name for d in datasets[d.dataset_status]]:
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
    endpoint = 'aggregate/all_runs'
    item_id = 'run_id'
    expected_bcl_subdirs = ('SampleSheet.csv', 'RunInfo.xml', 'Data')

    def __init__(self, config):
        super().__init__(config)
        self.use_int_dir = 'intermediate_dir' in config

    def get_dataset(self, name, most_recent_proc=None):
        dataset_path = os.path.join(self.input_dir, name)
        if os.path.exists(dataset_path):
            return RunDataset(
                name,
                os.path.join(self.input_dir, name),
                use_int_dir=self.use_int_dir,
                most_recent_proc=most_recent_proc
            )

    def _get_dataset_records_for_status(self, status):
        if status in (DATASET_NEW, DATASET_READY):
            self.debug('Scanning disk for datasets with status %s', status)
            datasets = []
            rest_api_datasets = rest_communication.get_documents(self.endpoint)
            for d in self._datasets_on_disk():
                if d not in rest_api_datasets:
                    datasets.append({self.item_id: d})
            return datasets
        else:
            return super()._get_dataset_records_for_status(status)

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
    endpoint = 'aggregate/samples'
    item_id = 'sample_id'

    def get_dataset(self, name, most_recent_proc=None, data_threshold=None):
        return SampleDataset(name, most_recent_proc, data_threshold)

    def _get_datasets_for_status(self, status):
        self.debug('Creating Datasets for status %s', status)
        datasets = {}
        for r in self._get_dataset_records_for_status(status):
            datasets[r[self.item_id]] = {'record': r}

        if datasets:
            self.debug('Querying Lims for %s samples', len(datasets))
            for sample in get_list_of_samples(list(datasets)):
                req_yield = sample.udf.get('Yield for Quoted Coverage (Gb)') * 1000000000
                datasets[sanitize_user_id(sample.name)]['threshold'] = req_yield
            self.debug('Lims query complete')

        return [
            self.get_dataset(k, v['record'].get('most_recent_proc'), v.get('threshold'))
            for k, v in datasets.items()
        ]
