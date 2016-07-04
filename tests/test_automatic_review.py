from bin import automatic_review

class TestAutomaticReview():


    @property
    def run_elements_by_lane_pass(self):
        run_elements_by_lane_pass = [{'lane_pc_optical_dups': None, 'sample_ids': ['LP1251760_G_09'], 'pc_q30': 83.86778972069591, 'lane_number': 7, 'cv': 0.0, 'stdev_pf': 0.0, 'useable': ['yes'], 'pc_pass_filter': 75.71414257913844, 'run_id': '160415_E00375_0091_BHFHWFCCXX', 'yield_q30_in_gb': 118.340301375, 'passing_filter_reads': 470344661, 'pc_q30_r2': 75.66213897769747, 'clean_yield_in_gb': 138.934823074, '_id': {'run_id': '160415_E00375_0091_BHFHWFCCXX', 'lane': 7}, 'avg_pf': 470344661.0, 'pc_q30_r1': 92.07344046369434, 'yield_in_gb': 141.1033983, 'clean_yield_q30_in_gb': 117.290061906},
                                {'lane_pc_optical_dups': None, 'sample_ids': ['LP1251760_A_10'], 'pc_q30': 83.01783234937763, 'lane_number': 1, 'cv': 0.0, 'stdev_pf': 0.0, 'useable': ['yes'], 'pc_pass_filter': 75.2795680870508, 'run_id': '160415_E00375_0091_BHFHWFCCXX', 'yield_q30_in_gb': 116.468631593, 'passing_filter_reads': 467645036, 'pc_q30_r2': 76.1873267769132, 'clean_yield_in_gb': 130.418415308, '_id': {'run_id': '160415_E00375_0091_BHFHWFCCXX', 'lane': 1}, 'avg_pf': 467645036.0, 'pc_q30_r1': 89.84833792184207, 'yield_in_gb': 140.2935108, 'clean_yield_q30_in_gb': 112.168373129},
                                {'lane_pc_optical_dups': None, 'sample_ids': ['LP1251760_B_10'], 'pc_q30': 81.42653323117273, 'lane_number': 2, 'cv': 0.0, 'stdev_pf': 0.0, 'useable': ['yes'], 'pc_pass_filter': 69.9352157748938, 'run_id': '160415_E00375_0091_BHFHWFCCXX', 'yield_q30_in_gb': 106.126130324, 'passing_filter_reads': 434445326, 'pc_q30_r2': 73.08344812376536, 'clean_yield_in_gb': 125.74650586, '_id': {'run_id': '160415_E00375_0091_BHFHWFCCXX', 'lane': 2}, 'avg_pf': 434445326.0, 'pc_q30_r1': 89.76961833858009, 'yield_in_gb': 130.3335978, 'clean_yield_q30_in_gb': 103.960923664},
                                {'lane_pc_optical_dups': None, 'sample_ids': ['LP1251760_D_10'], 'pc_q30': 82.29355551532952, 'lane_number': 4, 'cv': 0.0, 'stdev_pf': 0.0, 'useable': ['yes'], 'pc_pass_filter': 75.54566040081602, 'run_id': '160415_E00375_0091_BHFHWFCCXX', 'yield_q30_in_gb': 115.860610702, 'passing_filter_reads': 469298031, 'pc_q30_r2': 73.74754101621194, 'clean_yield_in_gb': 135.156955939, '_id': {'run_id': '160415_E00375_0091_BHFHWFCCXX', 'lane': 4}, 'avg_pf': 469298031.0, 'pc_q30_r1': 90.8395700144471, 'yield_in_gb': 140.7894093, 'clean_yield_q30_in_gb': 113.129957186},
                                {'lane_pc_optical_dups': None, 'sample_ids': ['LP1251760_C_10'], 'pc_q30': 82.29394041420912, 'lane_number': 3, 'cv': 0.0, 'stdev_pf': 0.0, 'useable': ['yes'], 'pc_pass_filter': 73.73193831383928, 'run_id': '160415_E00375_0091_BHFHWFCCXX', 'yield_q30_in_gb': 113.079524503, 'passing_filter_reads': 458030988, 'pc_q30_r2': 73.03492914035473, 'clean_yield_in_gb': 136.318772411, '_id': {'run_id': '160415_E00375_0091_BHFHWFCCXX', 'lane': 3}, 'avg_pf': 458030988.0, 'pc_q30_r1': 91.55295168806352, 'yield_in_gb': 137.4092964, 'clean_yield_q30_in_gb': 112.553373245},
                                {'lane_pc_optical_dups': None, 'sample_ids': ['LP1251760_F_09'], 'pc_q30': 83.28323516281932, 'lane_number': 6, 'cv': 0.0, 'stdev_pf': 0.0, 'useable': ['yes'], 'pc_pass_filter': 77.19009881059692, 'run_id': '160415_E00375_0091_BHFHWFCCXX', 'yield_q30_in_gb': 119.806298008, 'passing_filter_reads': 479513465, 'pc_q30_r2': 74.89543134866227, 'clean_yield_in_gb': 141.365054899, '_id': {'run_id': '160415_E00375_0091_BHFHWFCCXX', 'lane': 6}, 'avg_pf': 479513465.0, 'pc_q30_r1': 91.67103897697638, 'yield_in_gb': 143.8540395, 'clean_yield_q30_in_gb': 118.602469666},
                                {'lane_pc_optical_dups': None, 'sample_ids': ['LP1251760_H_09'], 'pc_q30': 85.08556857877406, 'lane_number': 8, 'cv': 0.0, 'stdev_pf': 0.0, 'useable': ['yes'], 'pc_pass_filter': 81.57124989188861, 'run_id': '160415_E00375_0091_BHFHWFCCXX', 'yield_q30_in_gb': 129.346144221, 'passing_filter_reads': 506729662, 'pc_q30_r2': 77.34754585835421, 'clean_yield_in_gb': 149.776224966, '_id': {'run_id': '160415_E00375_0091_BHFHWFCCXX', 'lane': 8}, 'avg_pf': 506729662.0, 'pc_q30_r1': 92.8235912991939, 'yield_in_gb': 152.0188986, 'clean_yield_q30_in_gb': 128.248403892},
                                {'lane_pc_optical_dups': None, 'sample_ids': ['LP1251760_E_09'], 'pc_q30': 83.27214896150255, 'lane_number': 5, 'cv': 0.0, 'stdev_pf': 0.0, 'useable': ['yes'], 'pc_pass_filter': 77.09395130193938, 'run_id': '160415_E00375_0091_BHFHWFCCXX', 'yield_q30_in_gb': 119.641139942, 'passing_filter_reads': 478916186, 'pc_q30_r2': 74.60413781602, 'clean_yield_in_gb': 141.870556438, '_id': {'run_id': '160415_E00375_0091_BHFHWFCCXX', 'lane': 5}, 'avg_pf': 478916186.0, 'pc_q30_r1': 91.94016010698512, 'yield_in_gb': 143.6748558, 'clean_yield_q30_in_gb': 118.760712074}]
        return run_elements_by_lane_pass

    @property
    def run_elements_by_lane_fail(self):
        run_elements_by_lane_fail = [{'useable': ['yes'], 'pc_q30': 71.61387317141885, 'stdev_pf': 0.0, 'yield_q30_in_gb': 52.673601996, 'cv': 0.0, 'clean_yield_q30_in_gb': 51.991327149, 'pc_pass_filter': 39.46711583571436, '_id': {'run_id': '160408_E00306_0095_AHFNCCCCXX', 'lane': 4}, 'yield_in_gb': 73.5522318, 'lane_pc_optical_dups': None, 'run_id': '160408_E00306_0095_AHFNCCCCXX', 'sample_ids': ['LP1251759_D_10'], 'avg_pf': 245174106.0, 'pc_q30_r2': 62.276207322372514, 'passing_filter_reads': 245174106, 'clean_yield_in_gb': 71.996565441, 'pc_q30_r1': 80.95153902046518, 'lane_number': 4},
                                     {'useable': ['yes'], 'pc_q30': 74.97437097775223, 'stdev_pf': 0.0, 'yield_q30_in_gb': 56.129839263, 'cv': 0.0, 'clean_yield_q30_in_gb': 55.561428357, 'pc_pass_filter': 40.17172703339186, '_id': {'run_id': '160408_E00306_0095_AHFNCCCCXX', 'lane': 2}, 'yield_in_gb': 74.8653687, 'lane_pc_optical_dups': None, 'run_id': '160408_E00306_0095_AHFNCCCCXX', 'sample_ids': ['LP1251759_B_10'], 'avg_pf': 249551229.0, 'pc_q30_r2': 64.26036718096067, 'passing_filter_reads': 249551229, 'clean_yield_in_gb': 73.564756728, 'pc_q30_r1': 85.68837477454379, 'lane_number': 2},
                                     {'useable': ['yes'], 'pc_q30': 74.37337631023922, 'stdev_pf': 0.0, 'yield_q30_in_gb': 54.625068958, 'cv': 0.0, 'clean_yield_q30_in_gb': 54.052712677, 'pc_pass_filter': 39.41068976770899, '_id': {'run_id': '160408_E00306_0095_AHFNCCCCXX', 'lane': 3}, 'yield_in_gb': 73.4470743, 'lane_pc_optical_dups': None, 'run_id': '160408_E00306_0095_AHFNCCCCXX', 'sample_ids': ['LP1251759_C_10'], 'avg_pf': 244823581.0, 'pc_q30_r2': 64.77634653447319, 'passing_filter_reads': 244823581, 'clean_yield_in_gb': 72.132320368, 'pc_q30_r1': 83.97040608600526, 'lane_number': 3},
                                     {'useable': ['yes'], 'pc_q30': 70.71992351661193, 'stdev_pf': 0.0, 'yield_q30_in_gb': 50.565242334, 'cv': 0.0, 'clean_yield_q30_in_gb': 49.330205084, 'pc_pass_filter': 38.36629359413382, '_id': {'run_id': '160408_E00306_0095_AHFNCCCCXX', 'lane': 6}, 'yield_in_gb': 71.5007028, 'lane_pc_optical_dups': None, 'run_id': '160408_E00306_0095_AHFNCCCCXX', 'sample_ids': ['LP1251759_F_10'], 'avg_pf': 238335676.0, 'pc_q30_r2': 61.80904496228252, 'passing_filter_reads': 238335676, 'clean_yield_in_gb': 68.793101585, 'pc_q30_r1': 79.63080207094131, 'lane_number': 6},
                                     {'useable': ['yes'], 'pc_q30': 63.71576692553782, 'stdev_pf': 0.0, 'yield_q30_in_gb': 21.363528373, 'cv': 0.0, 'clean_yield_q30_in_gb': 20.153848539, 'pc_pass_filter': 17.991426631034592, '_id': {'run_id': '160408_E00306_0095_AHFNCCCCXX', 'lane': 7}, 'yield_in_gb': 33.529422, 'lane_pc_optical_dups': None, 'run_id': '160408_E00306_0095_AHFNCCCCXX', 'sample_ids': ['LP1251759_G_10'], 'avg_pf': 111764740.0, 'pc_q30_r2': 51.25918433070513, 'passing_filter_reads': 111764740, 'clean_yield_in_gb': 30.603833777, 'pc_q30_r1': 76.1723495203705, 'lane_number': 7},
                                     {'useable': ['yes'], 'pc_q30': 70.56196810571603, 'stdev_pf': 0.0, 'yield_q30_in_gb': 48.627953338, 'cv': 0.0, 'clean_yield_q30_in_gb': 48.237385209, 'pc_pass_filter': 36.97897261025134, '_id': {'run_id': '160408_E00306_0095_AHFNCCCCXX', 'lane': 5}, 'yield_in_gb': 68.9152452, 'lane_pc_optical_dups': None, 'run_id': '160408_E00306_0095_AHFNCCCCXX', 'sample_ids': ['LP1251759_E_10'], 'avg_pf': 229717484.0, 'pc_q30_r2': 61.85159254718867, 'passing_filter_reads': 229717484, 'clean_yield_in_gb': 67.971434371, 'pc_q30_r1': 79.2723436642434, 'lane_number': 5},
                                     {'useable': ['yes'], 'pc_q30': 73.65629585810404, 'stdev_pf': 0.0, 'yield_q30_in_gb': 53.338724965, 'cv': 0.0, 'clean_yield_q30_in_gb': 50.898684919, 'pc_pass_filter': 38.85727049077345, '_id': {'run_id': '160408_E00306_0095_AHFNCCCCXX', 'lane': 8}, 'yield_in_gb': 72.4157037, 'lane_pc_optical_dups': None, 'run_id': '160408_E00306_0095_AHFNCCCCXX', 'sample_ids': ['LP1251759_H_10'], 'avg_pf': 241385679.0, 'pc_q30_r2': 60.31530744622178, 'passing_filter_reads': 241385679, 'clean_yield_in_gb': 66.770943534, 'pc_q30_r1': 86.99728426998632, 'lane_number': 8},
                                     {'useable': ['yes'], 'pc_q30': 77.14766136658056, 'stdev_pf': 0.0, 'yield_q30_in_gb': 57.325108732, 'cv': 0.0, 'clean_yield_q30_in_gb': 56.560546277, 'pc_pass_filter': 39.87141575627727, '_id': {'run_id': '160408_E00306_0095_AHFNCCCCXX', 'lane': 1}, 'yield_in_gb': 74.3056986, 'lane_pc_optical_dups': None, 'run_id': '160408_E00306_0095_AHFNCCCCXX', 'sample_ids': ['LP1251759_A_10'], 'avg_pf': 247685662.0, 'pc_q30_r2': 66.4934169100188, 'passing_filter_reads': 247685662, 'clean_yield_in_gb': 72.621861066, 'pc_q30_r1': 87.80190582314235, 'lane_number': 1}]
        return run_elements_by_lane_fail

    @property
    def samples_fail(self):
        samples_fail = {'properly_mapped_reads': 760174685,
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
          'duplicate_reads': 130972728,
          'run_ids': ['151030_E00365_0015_BHFCJ3CCXX'],
          'median_coverage': 30.156,
          'pc_properly_mapped_reads': 94.03266627960335,
          'clean_pc_q30_r1': 88.07789925713581,
          'clean_yield_q30': 89.332777716}
        return samples_fail


    @property
    def samples_pass(self):
        samples_pass = {'mapped_reads': 796025840,
         'library_id': 'LP0328602-DTP_B12',
         'sample_id': 'LP1251551__B_12',
         'genotype_validation': {'mismatching_snps': 1, 'no_call_seq': 0, 'no_call_chip': 3, 'matching_snps': 28},
         'proc_status': 'finished',
         'run_ids': ['151030_E00365_0014_AHFCNTCCXX'],
         'provided_gender': 'female',
         'clean_yield_in_gb': 120.173780067,
         'genotype_match': 'Match',
         'user_sample_id': 'LBC360090',
         'useable': 'yes',
         'pc_properly_mapped_reads': 92.71451664140578,
         'properly_mapped_reads': 771788247,
         'gender_match': 'female',
         '_id': 'LP1251551__B_12',
         'clean_yield_q30': 95.822039802,
         'called_gender': 'female',
         'clean_pc_q30_r1': 91.19261385088816,
         'median_coverage': 30.348,
         'duplicate_reads': 122371572,
         'pc_pass_filter': 100.0,
         'pc_mapped_reads': 95.62616595490734,
         'project_id': 'X0001ML',
         'clean_pc_q30_r2': 68.60350926684019,
         'delivered': 'yes',
         'pc_duplicate_reads': 20
         }

        return samples_pass










    def test_run_review(self):
        for lane in self.run_elements_by_lane_pass:
            run_review = automatic_review.metrics(lane, 'run')
            run_review = run_review[0]
            assert run_review == 'pass'
        for lane in self.run_elements_by_lane_fail:
            run_review = automatic_review.metrics(lane, 'run')
            print(run_review)
            assert run_review[0] == 'fail'



    def test_sample_review(self):
        sample_review = automatic_review.metrics(self.samples_fail, 'sample')
        review = sample_review[0]
        failed_fields = sample_review[1]
        assert review == 'fail'
        assert set(failed_fields) == set(['clean_yield_q30', 'clean_yield_in_gb'])


        sample_review = automatic_review.metrics(self.samples_pass, 'sample')
        review = sample_review[0]
        assert review == 'pass'

