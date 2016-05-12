import os
import sys
import argparse
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
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

    def lanes_metrics(self):
        lane_metrics = {}
        for r in self.run_element_by_lane:
            metrics_per_lane = {}
            lane = (r.get('lane_number'))
            metrics_per_lane['yield_in_gb'] = r.get('yield_in_gb')
            metrics_per_lane['pc_q30'] = r.get('pc_q30')
            metrics_per_lane['pc_pass_filter'] = r.get('pc_pass_filter')
            lane_metrics[lane] = metrics_per_lane
        return lane_metrics

    def run_review(self):
        lanes_metrics = self.lanes_metrics()
        run_review = {}
        for lane in lanes_metrics:
            lane_review = {}
            pc_pass_filter = (lanes_metrics[lane]).get('pc_pass_filter')
            yield_in_gb = (lanes_metrics[lane]).get('yield_in_gb')
            pc_q30 = (lanes_metrics[lane]).get('pc_q30')

            if float(yield_in_gb) > self.run_config.get('yield_in_gb'):
                lane_review['yield_in_gb'] = 'pass'
            else:
                lane_review['yield_in_gb'] = 'fail'

            if float(pc_q30) > self.run_config.get('pc_q30'):
                lane_review['pc_q30'] = 'pass'
            else:
                lane_review['pc_q30'] = 'fail'

            if float(pc_pass_filter) > self.run_config.get('pc_pass_filter'):
                lane_review['pc_pass_filter'] = 'pass'
            else:
                lane_review['pc_pass_filter'] = 'fail'

            if all(value == 'pass' for value in lane_review.values()):
                run_review[lane] = ['pass']
            else:
                failed_metrics = ([v for v in lane_review if lane_review[v] == 'fail'])
                run_review[lane] = ['fail', failed_metrics]
        return run_review

class SampleReview():

    def __init__(self, sample):
        self.sample = sample
        self.sample_cfg = cfg.get('sample')
    def sample_metrics(self):
        for r in self.sample:
            sample_id = (r.get('sample_id'))
            clean_yield_in_gb = (r.get('clean_yield_in_gb'))
            clean_yield_q30 = (r.get('clean_yield_q30'))
            pc_mapped_reads = (r.get('pc_mapped_reads'))
            median_coverage = (r.get('median_coverage'))
            pc_duplicate_reads = (r.get('pc_duplicate_reads'))
            genotype_validation = r.get('genotype_validation')
            called_gender = (r.get('called_gender'))
            provided_gender = (r.get('provided_gender'))
            return sample_id, \
                   clean_yield_in_gb, \
                   clean_yield_q30, \
                   pc_mapped_reads, \
                   median_coverage, \
                   pc_duplicate_reads, \
                   genotype_validation, \
                   called_gender, \
                   provided_gender

    def sample_review(self):
        sample_id,\
        clean_yield_in_gb, \
        clean_yield_q30, \
        pc_mapped_reads, \
        median_coverage, \
        pc_duplicate_reads, \
        genotype_validation, \
        called_gender, \
        provided_gender = self.sample_metrics()
        sample_metrics = {}
        sample_review = {}
        if clean_yield_in_gb and float(clean_yield_in_gb) > self.sample_cfg.get('clean_yield_in_gb'):
            sample_metrics['clean_yield_in_gb'] = 'pass'
        else:
            sample_metrics['clean_yield_in_gb'] = 'fail'

        if clean_yield_q30 and float(clean_yield_q30) > self.sample_cfg.get('clean_yield_q30'):
            sample_metrics['clean_yield_q30'] = 'pass'
        else:
            sample_metrics['clean_yield_q30'] = 'fail'

        if pc_mapped_reads and float(pc_mapped_reads) > self.sample_cfg.get('pc_mapped_reads'):
            sample_metrics['pc_mapped_reads'] = 'pass'
        else:
            sample_metrics['pc_mapped_reads'] = 'fail'

        if median_coverage and float(median_coverage) > self.sample_cfg.get('median_coverage'):
            sample_metrics['median_coverage'] = 'pass'
        else:
            sample_metrics['median_coverage'] = 'fail'

        if pc_duplicate_reads and float(pc_duplicate_reads) < self.sample_cfg.get('pc_duplicate_reads'):
            sample_metrics['pc_duplicate_reads'] = 'pass'
        else:
            sample_metrics['pc_duplicate_reads'] = 'fail'

        if provided_gender == 'unknown' or provided_gender == called_gender:
            sample_metrics['gender_validation'] = 'pass'
        else:
            sample_metrics['gender_validation'] = 'fail'

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
            sample_review[sample_id] = ['fail', failed_metrics]
        return sample_review

class AutomaticRunReview(RunReview):
    def __init__(self, run_name, run_element_by_lane):
        self.run_name = run_name
        self.run_element_by_lane = run_element_by_lane
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
    def __init__(self, sample_name, sample):
        self.sample = sample
        self.sample_name = sample_name
        super().__init__(sample)

    def patch_entry(self):
        sample_review = self.sample_review()
        review = ''
        if sample_review:
            if sample_review == 'pass':
                review = sample_review
            else:
                review = (list(sample_review.values())[0])[0]
            print(review)
            patch_review = {}
            patch_comments = {}

            if review == 'pass':
                patch_review['reviewed'] = 'pass'
                rest_communication.patch_entries('samples', payload=patch_review, update_lists=None, where={"sample_id": self.sample_name})

            else:
                patch_review['reviewed'] = 'fail'
                patch_comments[ELEMENT_REVIEW_COMMENTS] = (('sample ' + self.sample_name + ' failed due to ' + (' '.join((list(sample_review.values())[0])[1]))))
                rest_communication.patch_entries('samples', payload=patch_review, update_lists=None, where={"sample_id": self.sample_name})
                rest_communication.patch_entries('samples', payload=patch_comments, update_lists=None, where={"sample_id": self.sample_name})

def get_reviewable_runs():
    runs = rest_communication.get_documents('aggregate/all_runs', depaginate=True, match={"proc_status":"finished","review_statuses":"not%20reviewed"})
    for run in runs:
        run_id = run.get('run_id')
        run_elements_by_lane = rest_communication.get_documents('aggregate/run_elements_by_lane', match={"run_id":run_id})
        r = AutomaticRunReview(run_id, run_elements_by_lane)
        run_review = r.patch_entry()

def get_reviewable_samples():
    samples = rest_communication.get_documents('aggregate/samples', depaginate=True,  match={"proc_status":"finished","reviewed":"not%20reviewed"})
    for sp in samples:
        sample_id = sp.get('sample_id')
        sample = rest_communication.get_documents('aggregate/samples', match={"sample_id":sample_id})
        s = AutomaticSampleReview(sample_id, sample)
        sample_review = s.patch_entry()


if __name__ == '__main__':
    sys.exit(main())

