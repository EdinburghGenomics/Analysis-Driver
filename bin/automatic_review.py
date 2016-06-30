import os
import sys
import argparse
import logging
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from analysis_driver.app_logging import logging_default as log_cfg
log_cfg.default_level = logging.DEBUG
log_cfg.add_handler(logging.StreamHandler(stream=sys.stdout), logging.DEBUG)
from analysis_driver.external_data import rest_communication
from analysis_driver.constants import ELEMENT_REVIEW_COMMENTS
from analysis_driver.config import review_thresholds_config as cfg

def main():
    args = _parse_args()
    if not 'sample' not in sys.argv and 'run' not in sys.argv:
        print('no mode specified - use either --run_review or --sample_review')
        return 1
    else:
        if args.run_review:
            get_reviewable_runs()
        elif args.sample_review:
            get_reviewable_samples()


def _parse_args():
    parser = argparse.ArgumentParser()
    group = parser.add_mutually_exclusive_group()
    group.add_argument('--run_review', action='store_true')
    group.add_argument('--sample_review', action='store_true')
    return parser.parse_args()


def query(content, parts, top_level=None, ret_default=None):
    if top_level is None:
        top_level = content
    item = None
    for p in parts:
        item = top_level.get(p)
        if item != None:
            top_level = item
        else:
            return ret_default

    return item


def morethan(a, b):
    if a is not None:
        return a >= b
    else:
        return False

def lessthan(a, b):

    if a is not None:
        return a <= b
    else:
        return False

def inlist(a, b):
    if a is not None:
        return a in b

def metrics(dataset, datatype):
    PassFailDict = {}
    return_failed_metrics = None

    for metric in (cfg.get(datatype)):
        metric_value = (query(dataset, metric.split(',')))
        metric_name = (metric.split(',')[-1])
        if ((cfg.get(datatype)).get(metric)['comparison']) == '>':
            if morethan(metric_value, ((cfg.get(datatype)).get(metric)['value'])):
                PassFailDict[metric_name] = 'pass'
            else:
                PassFailDict[metric_name] = 'fail'
        elif ((cfg.get(datatype)).get(metric)['comparison']) == '<':
            if lessthan(metric_value, ((cfg.get(datatype)).get(metric)['value'])):
                PassFailDict[metric_name] = 'pass'
            else:
                PassFailDict[metric_name] = 'fail'
        elif ((cfg.get(datatype)).get(metric)['comparison']) == 'inlist':
            values = (cfg.get(datatype).get(metric)['value'])
            values_list = []
            for v in values:
                value = (query(dataset, [v], ret_default=[v]))
                values_list.append(''.join(value))
            check_metric = (list(set(values_list)))
            if inlist(metric_value, check_metric):
                PassFailDict[metric_name] = 'pass'
            else:
                PassFailDict[metric_name] = 'fail'

    if all(value == 'pass' for value in PassFailDict.values()):
        return ('pass', None)
    else:
        failed_metrics = ([value for value in PassFailDict if PassFailDict[value] == 'fail'])
        return_failed_metrics = ('fail', failed_metrics)
    return return_failed_metrics


class AutomaticRunReview():
    def __init__(self, run_name, run_element_by_lane):
        self.run_name = run_name
        self.run_element_by_lane = run_element_by_lane

    def patch_entry(self):
        for lane in self.run_element_by_lane:
            lane_review, reasons = metrics(lane, 'run')
            if lane_review:
                patch_review = {}
                patch_comments = {}
                if lane_review == 'pass':
                    patch_review['reviewed'] = 'pass'
                    rest_communication.patch_entries('run_elements', payload=patch_review, update_lists=None, where={"run_id":self.run_name, "lane": lane})
                elif lane_review == 'fail':
                    review_comments = ('lane ' + str(lane) + ' failed due to ' + ' '.join(reasons))
                    patch_review['reviewed'] = 'fail'
                    patch_comments[ELEMENT_REVIEW_COMMENTS] = review_comments
                    rest_communication.patch_entries('run_elements', payload=patch_review, update_lists=None, where={"run_id":self.run_name, "lane": (lane)})
                    rest_communication.patch_entries('run_elements', payload=patch_comments, update_lists=None, where={"run_id":self.run_name, "lane": (lane)})



class AutomaticSampleReview():
    def __init__(self, sample, sample_name):
        self.sample = sample
        self.sample_name = sample_name

    def patch_entry(self):
        sample_review, reasons = metrics(self.sample, 'sample')
        if sample_review:
            patch_review = {}
            patch_comments = {}
            if sample_review == 'pass':
                patch_review['reviewed'] = 'pass'
                rest_communication.patch_entries('samples', payload=patch_review, update_lists=None, where={"sample_id": self.sample_name})
            elif sample_review == 'fail':
                patch_review['reviewed'] = 'fail'
                patch_comments[ELEMENT_REVIEW_COMMENTS] = 'sample ' + self.sample_name + ' failed due to ' + ', '.join(reasons)
                rest_communication.patch_entries('samples', payload=patch_review, update_lists=None, where={"sample_id": self.sample_name})
                rest_communication.patch_entries('samples', payload=patch_comments, update_lists=None, where={"sample_id": self.sample_name})

def get_reviewable_runs():
    runs = rest_communication.get_documents('aggregate/all_runs', paginate=False, match={"proc_status":"finished","review_statuses":"not%20reviewed"})
    if runs:
        for run in runs:
            run_id = run.get('run_id')
            run_elements_by_lane = rest_communication.get_documents('aggregate/run_elements_by_lane', match={"run_id":run_id})
            r = AutomaticRunReview(run_id, run_elements_by_lane)
            r.patch_entry()

def get_reviewable_samples():
    samples = rest_communication.get_documents('aggregate/samples', paginate=False,  match={"proc_status":"finished","reviewed":"not%20reviewed"})
    if samples:
        for sample in samples:

            sample_id = sample.get('sample_id')
            s = AutomaticSampleReview(sample, sample_id)
            s.patch_entry()


if __name__ == '__main__':
    sys.exit(main())