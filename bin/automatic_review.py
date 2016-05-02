import os
import sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from analysis_driver import rest_communication
from analysis_driver.constants import ELEMENT_REVIEW_COMMENTS



class RunReview():
    def __init__(self, run_name):
        self.run_name = run_name
        self.run = self.read_run()
        self.lane = self.read_lanes()

    def read_run(self):
        return rest_communication.get_documents(
            'aggregate/all_runs',
            match={"run_id":self.run_name}
        )

    def read_lanes(self):
        return rest_communication.get_documents(
            'aggregate/run_elements',
            match={"run_id":self.run_name}
        )

    def run_metrics(self):
        for r in self.run:
            proc_status = str(r.get("proc_status"))
            return proc_status

    def lanes_metrics(self):
        lane_metrics = {}
        for r in self.lane:
            metrics_per_lane = {}
            lane = str(r.get('lane'))
            metrics_per_lane['yield_in_gb'] = str(r.get('yield_in_gb'))
            metrics_per_lane['pc_q30'] = str(r.get('pc_q30'))
            metrics_per_lane['pc_pass_filter'] = str(r.get('pc_pass_filter'))
            lane_metrics[lane] = metrics_per_lane
        return lane_metrics

    def run_review(self):
        if self.run and self.lane:
            proc_status = self.run_metrics()
            lanes_metrics = self.lanes_metrics()
            run_review = {}
            for lane in lanes_metrics:
                lane_review = {}
                pc_pass_filter = (lanes_metrics[lane]).get('pc_pass_filter')
                yield_in_gb = (lanes_metrics[lane]).get('yield_in_gb')
                pc_q30 = (lanes_metrics[lane]).get('pc_q30')

                if proc_status == 'finished':
                    lane_review['proc_status'] = 'pass'
                else:
                    lane_review['proc_status'] = 'fail'

                if float(yield_in_gb) > 960:
                    lane_review['yield_in_gb'] = 'pass'
                else:
                    lane_review['yield_in_gb'] = 'fail'

                if float(pc_q30) > 75:
                    lane_review['pc_q30'] = 'pass'
                else:
                    lane_review['pc_q30'] = 'fail'

                if float(pc_pass_filter) > 50:
                    lane_review['pc_pass_filter'] = 'pass'
                else:
                    lane_review['pc_pass_filter'] = 'fail'

                if all(value == 'pass' for value in lane_review.values()):
                    run_review[lane] = ['pass']
                else:
                    failed_metrics = ([v for v in lane_review if lane_review[v] == 'fail'])
                    run_review[lane] = ['fail', failed_metrics]
            return run_review
        return None


class SampleReview():

    def __init__(self, sample_name):
        self.sample_name = sample_name
        self.sample = self.read_sample()

    def read_sample(self):
        return rest_communication.get_documents(
            'aggregate/samples',
            match={"sample_id":self.sample_name}
        )

    def sample_metrics(self):
        for r in self.sample:
            proc_status = str(r.get('proc_status'))
            clean_yield_in_gb = str(r.get('clean_yield_in_gb'))
            clean_yield_q30 = str(r.get('clean_yield_q30'))
            pc_mapped_reads = str(r.get('pc_mapped_reads'))
            median_coverage = str(r.get('median_coverage'))
            pc_duplicate_reads = str(r.get('pc_duplicate_reads'))
            genotype_validation = r.get('genotype_validation')
            mismatching_snps = str(genotype_validation.get('mismatching_snps'))
            no_call_seq = str(genotype_validation.get('no_call_seq'))
            called_gender = str(r.get('called_gender'))
            provided_gender = str(r.get('provided_gender'))
            return proc_status, \
                   clean_yield_in_gb, \
                   clean_yield_q30, \
                   pc_mapped_reads, \
                   median_coverage, \
                   pc_duplicate_reads, \
                   mismatching_snps, \
                   no_call_seq, \
                   called_gender, \
                   provided_gender


    def sample_review(self):
        if self.sample:
            proc_status, \
            clean_yield_in_gb, \
            clean_yield_q30, \
            pc_mapped_reads, \
            median_coverage, \
            pc_duplicate_reads, \
            mismatching_snps, \
            no_call_seq, \
            called_gender, \
            provided_gender = self.sample_metrics()

            sample_metrics = {}
            sample_review = {}


            if proc_status != 'null' and proc_status == 'finished':
                sample_metrics['proc_status'] = 'pass'
            else:
                sample_metrics['proc_status'] = 'fail'

            if clean_yield_in_gb != 'null' and float(clean_yield_in_gb) > 120:
                sample_metrics['clean_yield_in_gb'] = 'pass'
            else:
                sample_metrics['clean_yield_in_gb'] = 'fail'

            if clean_yield_q30 != 'null' and float(clean_yield_q30) > 95:
                sample_metrics['clean_yield_q30'] = 'pass'
            else:
                sample_metrics['clean_yield_q30'] = 'fail'

            if pc_mapped_reads != 'null' and float(pc_mapped_reads) > 90:
                sample_metrics['pc_mapped_reads'] = 'pass'
            else:
                sample_metrics['pc_mapped_reads'] = 'fail'

            if median_coverage != 'null' and float(median_coverage) > 30:
                sample_metrics['median_coverage'] = 'pass'
            else:
                sample_metrics['median_coverage'] = 'fail'

            if pc_duplicate_reads != 'null' and float(pc_duplicate_reads) < 35:
                sample_metrics['pc_duplicate_reads'] = 'pass'
            else:
                sample_metrics['pc_duplicate_reads'] = 'fail'

            if provided_gender == 'unknown' or provided_gender == called_gender:
                sample_metrics['gender_validation'] = 'pass'
            else:
                sample_metrics['gender_validation'] = 'fail'

            if mismatching_snps != 'null' and int(mismatching_snps) < 5:
                sample_metrics['mismatching_snps'] = 'pass'
            else:
                sample_metrics['mismatching_snps'] = 'fail'

            if no_call_seq != 'null' and int(no_call_seq) < 10:
                sample_metrics['no_call_seq'] = 'pass'
            else:
                sample_metrics['no_call_seq'] = 'fail'



            if all(value == 'pass' for value in sample_metrics.values()):
                return 'pass'
            else:
                failed_metrics = ([value for value in sample_metrics if sample_metrics[value] == 'fail'])
                sample_review[self.sample_name] = ['fail', failed_metrics]
            return sample_review
        return None


class AutomaticRunReview(RunReview):
    def __init__(self, run_name):
        super().__init__(run_name)

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
                    review_comments = ('lane ' + lane + ' failed due to ' + ' '.join(run_review[lane][1]))

                    patch_review['reviewed'] = 'fail'
                    patch_comments[ELEMENT_REVIEW_COMMENTS] = review_comments
                    rest_communication.patch_entries('run_elements', payload=patch_review, update_lists=None, where={"run_id":self.run_name, "lane": lane})
                    rest_communication.patch_entries('run_elements', payload=patch_comments, update_lists=None, where={"run_id":self.run_name, "lane": lane})



class AutomaticSampleReview(SampleReview):
    def __init__(self, sample_name):
        super().__init__(sample_name)

    def patch_entry(self):
        sample_review = self.sample_review()
        if sample_review:

            review = (list(sample_review.values())[0])[0]

            patch_review = {}
            patch_comments = {}

            if review == 'pass':
                patch_review['reviewed'] = 'pass'
                rest_communication.patch_entries('samples', payload=patch_review, update_lists=None, where={"sample_id": self.sample_name})

            else:
                patch_review['reviewed'] = 'fail'
                patch_comments[ELEMENT_REVIEW_COMMENTS] = ('sample '+ self.sample_name + ' failed due to ' + ' '.join(list(sample_review.values())[0])[0])
                rest_communication.patch_entries('samples', payload=patch_review, update_lists=None, where={"sample_id": self.sample_name})
                rest_communication.patch_entries('samples', payload=patch_comments, update_lists=None, where={"sample_id": self.sample_name})






s = AutomaticSampleReview('LP1251551__C_04')
sample_review = s.patch_entry()

