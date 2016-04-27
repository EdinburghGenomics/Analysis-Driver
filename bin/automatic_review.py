import os
import sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from analysis_driver import rest_communication



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
            'aggregate/run_elements_by_lane',
            match={"run_id":self.run_name}
        )

    def run_metrics(self):
        for r in self.run:
            proc_status = str(r.get("proc_status"))
            #clean_yield_in_gb = str(r.get("clean_yield_in_gb"))
            #pc_q30 = str(r.get("pc_q30"))
            return proc_status

    def lanes_metrics(self):
        for r in self.lane:
            lane = str(r.get('lane_number'))
            yield_in_gb = str(r.get('yield_in_gb'))
            pc_q30 = str(r.get('pc_q30'))
            pc_pass_filter = str(r.get('pc_pass_filter'))

            return lane, yield_in_gb, pc_q30, pc_pass_filter

    def run_review(self):
        if self.run and self.lane:
            proc_status = self.run_metrics()
            lane, \
            yield_in_gb, \
            pc_q30, \
            pc_pass_filter = self.lanes_metrics()

            if proc_status == 'finished':
                if float(yield_in_gb) > 960:
                    if float(pc_q30) > 75:
                        if all(float(v) > 50 for v in pcPassFilterByLane.values()):
                            return 'review pass'
            else:
                return 'review fail'
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
            pc_properly_mapped_reads = str(r.get('pc_properly_mapped_reads'))
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
                   pc_properly_mapped_reads, \
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
            pc_properly_mapped_reads, \
            median_coverage, \
            pc_duplicate_reads, \
            mismatching_snps, \
            no_call_seq, \
            called_gender, \
            provided_gender = self.sample_metrics()

            if proc_status == 'finished':
                if float(clean_yield_in_gb) > 120:
                    if float(clean_yield_q30) > 95:
                        if float(pc_properly_mapped_reads) > 90:
                            if float(median_coverage) > 30:
                                if pc_duplicate_reads < 35:
                                    if provided_gender == 'unknown' or provided_gender == called_gender:
                                        if int(mismatching_snps) < 5:
                                            if no_call_seq < 10:
                                                return 'review pass'
            else:
                return 'review fail'
        return None







r = RunReview("150723_E00306_0025_BHCHK3CCXX")
runReview = r.run_review()

s = SampleReview("LP1")
sampleReview = s.sample_review()


print(runReview)
print(sampleReview)

#rest_communication.patch_entries(end_point, payload=patch, update_lists=None, where=f)