import os
from unittest.mock import patch
from bin import automatic_review


class TestAutomaticReview():

    @property
    def run_elements_by_lane(self):
        run_elements_by_lane = [{'cv': 0.39289755963001993,
                                 'avg_pf': 70076060.57142857,
                                 'pc_pass_filter': 78.96388535900994,
                                 'pc_q30': 88.3062505602416,
                                 'stdev_pf': 27532713.186999746,
                                 'lane_number': 8,
                                 'yield_q30_in_gb': 129.951237425,
                                 'pc_q30_r1': 95.71726501950188,
                                 'pc_q30_r2': 80.89523610098131,
                                 'passing_filter_reads': 490532424,
                                 'yield_in_gb': 147.1597272,
                                 '_id': 8},
                                {'cv': 0.3985219296212109,
                                 'avg_pf': 69605637.85714285,
                                 'pc_pass_filter': 78.43379840808512,
                                 'pc_q30': 88.30995685116216,
                                 'stdev_pf': 27739373.11134378,
                                 'lane_number': 6,
                                 'yield_q30_in_gb': 129.084288391,
                                 'pc_q30_r1': 95.60978256417168,
                                 'pc_q30_r2': 81.01013113815264,
                                 'passing_filter_reads': 487239465,
                                 'yield_in_gb': 146.1718395,
                                 '_id': 6},
                                {'cv': 0.4118516395030433,
                                 'avg_pf': 67664981.0,
                                 'pc_pass_filter': 76.24700588095091,
                                 'pc_q30': 88.64660997490958,
                                 'stdev_pf': 27867933.361792278,
                                 'lane_number': 5,
                                 'yield_q30_in_gb': 125.963694773,
                                 'pc_q30_r1': 95.32747929446836,
                                 'pc_q30_r2': 81.96574065535079,
                                 'passing_filter_reads': 473654867,
                                 'yield_in_gb': 142.0964601,
                                 '_id': 5},
                                {'cv': 0.4033879132046018,
                                 'avg_pf': 68325216.85714285,
                                 'pc_pass_filter': 76.99098018698648,
                                 'pc_q30': 88.07618521913997,
                                 'stdev_pf': 27561566.64725474,
                                 'lane_number': 4,
                                 'yield_q30_in_gb': 126.374313556,
                                 'pc_q30_r1': 95.25712949316627,
                                 'pc_q30_r2': 80.89524094511368,
                                 'passing_filter_reads': 478276518,
                                 'yield_in_gb': 143.4829554,
                                 '_id': 4},
                                {'cv': 0.41549048201463856,
                                 'avg_pf': 67341506.71428572,
                                 'pc_pass_filter': 75.88250499141111,
                                 'pc_q30': 88.23926437017273,
                                 'stdev_pf': 27979755.08431059,
                                 'lane_number': 3,
                                 'yield_q30_in_gb': 124.785465295,
                                 'pc_q30_r1': 94.84886269473706,
                                 'pc_q30_r2': 81.62966604560839,
                                 'passing_filter_reads': 471390547,
                                 'yield_in_gb': 141.4171641,
                                 '_id': 3},
                                {'cv': 0.4190527974338691,
                                 'avg_pf': 66747671.0,
                                 'pc_pass_filter': 75.21335243228363,
                                 'pc_q30': 88.39428342072968,
                                 'stdev_pf': 27970798.25474554,
                                 'lane_number': 2,
                                 'yield_q30_in_gb': 123.902363509,
                                 'pc_q30_r1': 95.18075267731956,
                                 'pc_q30_r2': 81.6078141641398,
                                 'passing_filter_reads': 467233697,
                                 'yield_in_gb': 140.1701091,
                                 '_id': 2},
                                {'cv': 0.39546625347962017,
                                 'avg_pf': 69094749.14285715,
                                 'pc_pass_filter': 77.85811310932395,
                                 'pc_q30': 89.09804778825271,
                                 'stdev_pf': 27324641.578639913,
                                 'lane_number': 7,
                                 'yield_q30_in_gb': 129.280352482,
                                 'pc_q30_r1': 95.70192734072359,
                                 'pc_q30_r2': 82.49416823578184,
                                 'passing_filter_reads': 483663244,
                                 'yield_in_gb': 145.0989732,
                                 '_id': 7},
                                {'cv': 0.4155146587647161,
                                 'avg_pf': 63743309.571428575,
                                 'pc_pass_filter': 71.8279444985581,
                                 'pc_q30': 88.72861682758966,
                                 'stdev_pf': 26486279.525105808,
                                 'lane_number': 1,
                                 'yield_q30_in_gb': 118.772969496,
                                 'pc_q30_r1': 94.80184001473033,
                                 'pc_q30_r2': 82.65539364044899,
                                 'passing_filter_reads': 446203167,
                                 'yield_in_gb': 133.8609501,
                                 '_id': 1}]
        return run_elements_by_lane

    @property
    def samples(self):
        samples = [{'properly_mapped_reads': 760174685,
          'delivered': 'yes',
          'pc_mapped_reads': 97.18554753953777,
          'review_comments': 'sample LP1251551__C_04 failed due to clean_yield_in_gb clean_yield_q30',
          'proc_status': 'finished',
          'reviewed': 'fail',
          'clean_pc_q30': 78.14874745012538,
          'bam_file_reads': 808415538,
          'user_sample_id': 'LBC360021',
          'provided_gender': 'male',
          'clean_yield_in_gb': 114.3112086,
          'sample_id': 'LP1251551__C_04',
          'project_id': 'X0001ML',
          '_id': 'LP1251551__C_04',
          'called_gender': 'male',
          'mapped_reads': 785663067,
          'useable': 'yes',
          'clean_pc_q30_r2': 68.22116840186138,
          'genotype_validation': {'mismatching_snps': 3, 'no_call_chip': 1, 'no_call_seq': 0, 'matching_snps': 28},
          'pc_duplicate_reads': 16.20116410973783,
          'pc_pass_filter': 100.0,
          'library_id': 'LP0328602-DTP_C04',
          'duplicate_reads': 130972728, 'run_ids': ['151030_E00365_0015_BHFCJ3CCXX'],
          'median_coverage': 30.156, 'pc_properly_mapped_reads': 94.03266627960335,
          'clean_pc_q30_r1': 88.07789925713581,
          'clean_yield_q30': 89.332777716}]
        return samples



    @property
    def assets_path(self):
        assets_path = os.path.join(os.path.dirname(__file__), 'assets')
        return assets_path


    def test_lanes_metrics(self):
        run_elements_by_lane = self.run_elements_by_lane
        r = automatic_review.RunReview(run_elements_by_lane)
        lanes_metrics = r.lanes_metrics()
        assert lanes_metrics == {'1': {'pc_pass_filter': '71.8279444985581', 'pc_q30': '88.72861682758966', 'yield_in_gb': '133.8609501'},
                                 '8': {'pc_pass_filter': '78.96388535900994', 'pc_q30': '88.3062505602416', 'yield_in_gb': '147.1597272'},
                                 '4': {'pc_pass_filter': '76.99098018698648', 'pc_q30': '88.07618521913997', 'yield_in_gb': '143.4829554'},
                                 '5': {'pc_pass_filter': '76.24700588095091', 'pc_q30': '88.64660997490958', 'yield_in_gb': '142.0964601'},
                                 '6': {'pc_pass_filter': '78.43379840808512', 'pc_q30': '88.30995685116216', 'yield_in_gb': '146.1718395'},
                                 '7': {'pc_pass_filter': '77.85811310932395', 'pc_q30': '89.09804778825271', 'yield_in_gb': '145.0989732'},
                                 '3': {'pc_pass_filter': '75.88250499141111', 'pc_q30': '88.23926437017273', 'yield_in_gb': '141.4171641'},
                                 '2': {'pc_pass_filter': '75.21335243228363', 'pc_q30': '88.39428342072968', 'yield_in_gb': '140.1701091'}}


    def test_run_review(self):
        run_elements_by_lane = self.run_elements_by_lane
        r = automatic_review.RunReview(run_elements_by_lane)
        run_review  = r.run_review()
        assert run_review == {'8': ['fail', ['yield_in_gb']],
                              '4': ['fail', ['yield_in_gb']],
                              '1': ['fail', ['yield_in_gb']],
                              '6': ['fail', ['yield_in_gb']],
                              '5': ['fail', ['yield_in_gb']],
                              '7': ['fail', ['yield_in_gb']],
                              '3': ['fail', ['yield_in_gb']],
                              '2': ['fail', ['yield_in_gb']]}



