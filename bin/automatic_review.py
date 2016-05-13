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



class RunReview():
    def __init__(self, run_element_by_lane):
        self.run_element_by_lane = run_element_by_lane
        self.run_config = cfg.get('run')


    def run_review(self):
        run_metrics = {}
        for lane in self.run_element_by_lane:
            lane_review = {}
            lane_id = lane['lane_number']

            if float(lane['yield_in_gb']) > self.run_config.get('yield_in_gb'):
                lane_review['yield_in_gb'] = 'pass'
            else:
                lane_review['yield_in_gb'] = 'fail'

            if float(lane['pc_q30']) > self.run_config.get('pc_q30'):
                lane_review['pc_q30'] = 'pass'
            else:
                lane_review['pc_q30'] = 'fail'

            if float(lane['pc_pass_filter']) > self.run_config.get('pc_pass_filter'):
                lane_review['pc_pass_filter'] = 'pass'
            else:
                lane_review['pc_pass_filter'] = 'fail'

            if all(value == 'pass' for value in lane_review.values()):
                run_metrics[lane_id] = ['pass']
            else:
                failed_metrics = ([v for v in lane_review if lane_review[v] == 'fail'])
                run_metrics[lane_id] = ['fail', failed_metrics]
        return run_metrics








class SampleReview():

    def __init__(self, sample, sample_name):
        self.s = sample
        self.sample_name = sample_name
        self.sample_cfg = cfg.get('sample')

    def sample_review(self):
        sample_metrics = {}
        sample_review = {}
        if self.s.get('clean_yield_in_gb') and float(self.s.get('clean_yield_in_gb')) > self.sample_cfg.get('clean_yield_in_gb'):
            sample_metrics['clean_yield_in_gb'] = 'pass'
        else:
            sample_metrics['clean_yield_in_gb'] = 'fail'

        if self.s.get('clean_yield_q30') and float(self.s.get('clean_yield_q30')) > self.sample_cfg.get('clean_yield_q30'):
            sample_metrics['clean_yield_q30'] = 'pass'
        else:
            sample_metrics['clean_yield_q30'] = 'fail'

        if self.s.get('pc_mapped_reads') and float(self.s.get('pc_mapped_reads')) > self.sample_cfg.get('pc_mapped_reads'):
            sample_metrics['pc_mapped_reads'] = 'pass'
        else:
            sample_metrics['pc_mapped_reads'] = 'fail'

        if self.s.get('median_coverage') and float(self.s.get('median_coverage')) > self.sample_cfg.get('median_coverage'):
            sample_metrics['median_coverage'] = 'pass'
        else:
            sample_metrics['median_coverage'] = 'fail'

        if self.s.get('pc_duplicate_reads') and float(self.s.get('pc_duplicate_reads')) < self.sample_cfg.get('pc_duplicate_reads'):
            sample_metrics['pc_duplicate_reads'] = 'pass'
        else:
            sample_metrics['pc_duplicate_reads'] = 'fail'

        if self.s.get('provided_gender') == 'unknown' or self.s.get('provided_gender') == self.s.get('called_gender'):
            sample_metrics['gender_validation'] = 'pass'
        else:
            sample_metrics['gender_validation'] = 'fail'

        genotype_validation = self.s.get('genotype_validation')
        if genotype_validation != None:
            mismatching_snps = str(genotype_validation.get('mismatching_snps'))
            if mismatching_snps != 'null' and int(mismatching_snps) < self.sample_cfg.get('mismatching_snps'):
                sample_metrics['mismatching_snps'] = 'pass'
            else:
                sample_metrics['mismatching_snps'] = 'fail'

            no_call_seq = str(genotype_validation.get('no_call_seq'))
            if no_call_seq != 'null' and int(no_call_seq) < self.sample_cfg.get('no_call_seq'):
                sample_metrics['no_call_seq'] = 'pass'
            else:
                sample_metrics['no_call_seq'] = 'fail'


        if all(value == 'pass' for value in sample_metrics.values()):
            return 'pass'
        else:
            failed_metrics = ([value for value in sample_metrics if sample_metrics[value] == 'fail'])
            sample_review[self.sample_name] = ['fail', failed_metrics]
        return sample_review


class AutomaticRunReview(RunReview):
    def __init__(self, run_name, run_element_by_lane):
        self.run_name = run_name
        super().__init__(run_element_by_lane)

    def patch_entry(self):
        run_review = self.run_review()
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

class AutomaticSampleReview(SampleReview):
    def __init__(self, sample, sample_name):
        super().__init__(sample, sample_name)

    def patch_entry(self):
        sample_review = self.sample_review()
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

