import os
import sys
import argparse
import logging
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from analysis_driver.app_logging import logging_default as log_cfg
log_cfg.default_level = logging.DEBUG
log_cfg.add_handler(logging.StreamHandler(stream=sys.stdout), logging.DEBUG)
from analysis_driver import rest_communication
from analysis_driver.constants import ELEMENT_REVIEW_COMMENTS
from analysis_driver.config import review_thresholds_config as cfg

def main():
    args = _parse_args()
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


def run_metrics(run_element_by_lane):
    metrics = {}
    run_config = cfg.get('run')
    for lane in run_element_by_lane:
        lane_review = {}
        lane_id = lane['lane_number']

        if float(lane['yield_in_gb']) > run_config.get('yield_in_gb'):
            lane_review['yield_in_gb'] = 'pass'
        else:
            lane_review['yield_in_gb'] = 'fail'

        if float(lane['pc_q30']) > run_config.get('pc_q30'):
            lane_review['pc_q30'] = 'pass'
        else:
            lane_review['pc_q30'] = 'fail'

        if float(lane['pc_pass_filter']) > run_config.get('pc_pass_filter'):
            lane_review['pc_pass_filter'] = 'pass'
        else:
            lane_review['pc_pass_filter'] = 'fail'

        if all(value == 'pass' for value in lane_review.values()):
            metrics[lane_id] = ['pass']
        else:
            failed_metrics = ([v for v in lane_review if lane_review[v] == 'fail'])
            metrics[lane_id] = ['fail', failed_metrics]
    return metrics


def sample_metrics(s, sample_name):
    metric_values = {}
    metrics = {}
    sample_cfg = cfg.get('sample')

    if s.get('clean_yield_in_gb') and float(s.get('clean_yield_in_gb')) > sample_cfg.get('clean_yield_in_gb'):
        metric_values['clean_yield_in_gb'] = 'pass'
    else:
        metric_values['clean_yield_in_gb'] = 'fail'

    if s.get('clean_yield_q30') and float(s.get('clean_yield_q30')) > sample_cfg.get('clean_yield_q30'):
        metric_values['clean_yield_q30'] = 'pass'
    else:
        metric_values['clean_yield_q30'] = 'fail'

    if s.get('pc_mapped_reads') and float(s.get('pc_mapped_reads')) > sample_cfg.get('pc_mapped_reads'):
        metric_values['pc_mapped_reads'] = 'pass'
    else:
        metric_values['pc_mapped_reads'] = 'fail'

    if s.get('median_coverage') and float(s.get('median_coverage')) > sample_cfg.get('median_coverage'):
        metric_values['median_coverage'] = 'pass'
    else:
        metric_values['median_coverage'] = 'fail'

    if s.get('pc_duplicate_reads') and float(s.get('pc_duplicate_reads')) < sample_cfg.get('pc_duplicate_reads'):
        metric_values['pc_duplicate_reads'] = 'pass'
    else:
        metric_values['pc_duplicate_reads'] = 'fail'

    if s.get('provided_gender') == 'unknown' or s.get('provided_gender') == s.get('called_gender'):
        metric_values['gender_validation'] = 'pass'
    else:
        metric_values['gender_validation'] = 'fail'

    genotype_validation = s.get('genotype_validation')
    if genotype_validation != None:
        mismatching_snps = str(genotype_validation.get('mismatching_snps'))
        if mismatching_snps != 'null' and int(mismatching_snps) < sample_cfg.get('mismatching_snps'):
            metric_values['mismatching_snps'] = 'pass'
        else:
            metric_values['mismatching_snps'] = 'fail'

        no_call_seq = str(genotype_validation.get('no_call_seq'))
        if no_call_seq != 'null' and int(no_call_seq) < sample_cfg.get('no_call_seq'):
            metric_values['no_call_seq'] = 'pass'
        else:
            metric_values['no_call_seq'] = 'fail'


    if all(value == 'pass' for value in metric_values.values()):
        return 'pass'
    else:
        failed_metrics = ([value for value in metric_values if metric_values[value] == 'fail'])
        metrics[sample_name] = ['fail', failed_metrics]
    return metrics


class AutomaticRunReview():
    def __init__(self, run_name, run_element_by_lane):
        self.run_name = run_name
        self.run_element_by_lane = run_element_by_lane

    def patch_entry(self):
        run_review = run_metrics(self.run_element_by_lane)
        if run_review:
            for lane in run_review:
                patch_review = {}
                patch_comments = {}

                review = run_review[lane][0]

                if review == 'pass':
                    patch_review['reviewed'] = 'pass'
                    rest_communication.patch_entries('run_elements', payload=patch_review, update_lists=None, where={"run_id":self.run_name, "lane": lane})
                elif review == 'fail':
                    review_comments = ('lane ' + str(lane) + ' failed due to ' + ' '.join(run_review[lane][1]))
                    patch_review['reviewed'] = 'fail'
                    patch_comments[ELEMENT_REVIEW_COMMENTS] = review_comments
                    rest_communication.patch_entries('run_elements', payload=patch_review, update_lists=None, where={"run_id":self.run_name, "lane": (lane)})
                    rest_communication.patch_entries('run_elements', payload=patch_comments, update_lists=None, where={"run_id":self.run_name, "lane": (lane)})

class AutomaticSampleReview():
    def __init__(self, sample, sample_name):
        self.sample = sample
        self.sample_name = sample_name

    def patch_entry(self):
        sample_review = sample_metrics(self.sample, self.sample_name)
        review = ''
        if sample_review:
            if sample_review == 'pass':
                review = sample_review
            else:
                review = list(sample_review.values())[0][0]
            patch_review = {}
            patch_comments = {}

            if review == 'pass':
                patch_review['reviewed'] = 'pass'
                rest_communication.patch_entries('samples', payload=patch_review, update_lists=None, where={"sample_id": self.sample_name})

            else:
                patch_review['reviewed'] = 'fail'
                patch_comments[ELEMENT_REVIEW_COMMENTS] = 'sample ' + self.sample_name + ' failed due to ' + ' '.join(list(sample_review.values())[0])[1]
                rest_communication.patch_entries('samples', payload=patch_review, update_lists=None, where={"sample_id": self.sample_name})
                rest_communication.patch_entries('samples', payload=patch_comments, update_lists=None, where={"sample_id": self.sample_name})

def get_reviewable_runs():
    runs = rest_communication.get_documents('aggregate/all_runs', depaginate=True, match={"proc_status":"finished","review_statuses":"not%20reviewed"})
    if runs:
        for run in runs:
            run_id = run.get('run_id')
            run_elements_by_lane = rest_communication.get_documents('aggregate/run_elements_by_lane', match={"run_id":run_id})
            r = AutomaticRunReview(run_id, run_elements_by_lane)
            run_review = r.patch_entry()

def get_reviewable_samples():
    samples = rest_communication.get_documents('aggregate/samples', depaginate=True,  match={"proc_status":"finished","reviewed":"not%20reviewed"})
    if samples:
        for sample in samples:
            sample_id = sample.get('sample_id')
            s = AutomaticSampleReview(sample, sample_id)
            sample_review = s.patch_entry()


if __name__ == '__main__':
    sys.exit(main())

