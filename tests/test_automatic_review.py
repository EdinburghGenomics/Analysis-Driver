from bin import automatic_review

class TestAutomaticReview():


    @property
    def run_elements_by_lane(self):
        run_elements_by_lane = [{"passing_filter_reads": 476925678,
        "pc_q30_r1": 89.81539119812291,
        "pc_q30_r2": 72.365607698173335,
        "pc_pass_filter": 76.773527538232798,
        "stdev_pf": 0.0,
        "yield_in_gb": 143.07770339999999,
        "yield_q30_in_gb": 116.022424286,
        "pc_q30": 81.090499448148108,
        "lane_number": 3,
        "_id": 3,
        "cv": 0.0,
        "avg_pf": 476925678.0},
        {"passing_filter_reads": 489534556,
        "pc_q30_r1": 86.894980333114631,
        "pc_q30_r2": 70.399391353011382,
        "pc_pass_filter": 78.803252686223715,
        "stdev_pf": 0.0,
        "yield_in_gb": 146.86036680000001,
        "yield_q30_in_gb": 115.501545607,
        "pc_q30": 78.647185843062999,
        "lane_number": 8,
        "_id": 8,
        "cv": 0.0,
        "avg_pf": 489534556.0},
        {"passing_filter_reads": 471087617,
        "pc_q30_r1": 89.505951924296355,
        "pc_q30_r2": 71.944623808695866,
        "pc_pass_filter": 75.833740570097731,
        "stdev_pf": 0.0,
        "yield_in_gb": 141.32628510000001,
        "yield_q30_in_gb": 114.086050478,
        "pc_q30": 80.725287866496103,
        "lane_number": 1,
        "_id": 1,
        "cv": 0.0,
        "avg_pf": 471087617.0},
        {"passing_filter_reads": 477352843,
        "pc_q30_r1": 84.748950282603985,
        "pc_q30_r2": 66.472150973097555,
        "pc_pass_filter": 76.842290797171586,
        "stdev_pf": 0.0,
        "yield_in_gb": 143.2058529,
        "yield_q30_in_gb": 108.278733909,
        "pc_q30": 75.61055062785077,
        "lane_number": 7,
        "_id": 7,
        "cv": 0.0,
        "avg_pf": 477352843.0},
        {"passing_filter_reads": 468438741,
        "pc_q30_r1": 85.036875527480476,
        "pc_q30_r2": 66.995624157111862,
        "pc_pass_filter": 75.407335442606637,
        "stdev_pf": 0.0,
        "yield_in_gb": 140.53162230000001,
        "yield_q30_in_gb": 106.82686911499999,
        "pc_q30": 76.016249842296162,
        "lane_number": 5,
        "_id": 5,
        "cv": 0.0,
        "avg_pf": 468438741.0},
        {"passing_filter_reads": 459498349,
        "pc_q30_r1": 86.713012366739392,
        "pc_q30_r2": 66.396725915838559,
        "pc_pass_filter": 73.968148032331371,
        "stdev_pf": 0.0,
        "yield_in_gb": 137.84950470000001,
        "yield_q30_in_gb": 105.530507935,
        "pc_q30": 76.554869141288989,
        "lane_number": 2,
        "_id": 2,
        "cv": 0.0,
        "avg_pf": 459498349.0},
        {"passing_filter_reads": 430231698,
        "pc_q30_r1": 86.973534282295788,
        "pc_q30_r2": 71.103524699691775,
        "pc_pass_filter": 69.25692332762938,
        "stdev_pf": 0.0,
        "yield_in_gb": 129.06950939999999,
        "yield_q30_in_gb": 102.014642251,
        "pc_q30": 79.038529490993795,
        "lane_number": 4,
        "_id": 4,
        "cv": 0.0,
        "avg_pf": 430231698.0},
        {"passing_filter_reads": 345154480,
        "pc_q30_r1": 86.184789577891806,
        "pc_q30_r2": 68.261867999897703,
        "pc_pass_filter": 55.561543858044104,
        "stdev_pf": 0.0,
        "yield_in_gb": 103.546344,
        "yield_q30_in_gb": 79.961933676000001,
        "pc_q30": 77.223328788894761,
        "lane_number": 6,
        "_id": 6,
        "cv": 0.0,
        "avg_pf": 345154480.0}]
        return run_elements_by_lane


    @property
    def failed_samples(self):
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
          'duplicate_reads': 130972728,
          'run_ids': ['151030_E00365_0015_BHFCJ3CCXX'],
          'median_coverage': 30.156,
          'pc_properly_mapped_reads': 94.03266627960335,
          'clean_pc_q30_r1': 88.07789925713581,
          'clean_yield_q30': 89.332777716}]
        return samples

    @property
    def passed_samples(self):
        samples = [{"library_id": "LP6002104-DTP_B06",
        "genotype_validation": {"matching_snps": 28, "no_call_chip": 2, "mismatching_snps":2, "no_call_seq": 0},
        "called_gender": "male",
        "reviewed": "pass",
        "provided_gender": "male",
        "sample_id": "LP1251553__B_06",
        "duplicate_reads": 268712798,
        "delivered": "yes",
        "user_sample_id": "LBC360165",
        "pc_duplicate_reads": 30.956808538365465,
        "pc_mapped_reads": 98.224465279810019,
        "clean_pc_q30": 79.1101906289298,
        "mapped_reads": 852612790,
        "median_coverage": 26.716000000000001,
        "project_id": "X0001ML",
        "clean_yield_q30": 127.51388991499999,
        "clean_yield_in_gb": 161.18516325300001,
        "run_ids": ["160128_E00306_0061_BHFGCCCCXX", "160121_E00328_0031_AHFHLLCCXX"],
        "clean_pc_q30_r2": 68.785578424863559,
        "clean_pc_q30_r1": 89.438872446261612,
        "useable": "yes",
        "bam_file_reads": 868024873,
        "pc_pass_filter": 100.0,
        "review_comments": None,
        "pc_properly_mapped_reads": 95.066261655384608,
        "properly_mapped_reads": 825198797,
        "_id": "LP1251553__B_06",
        "proc_status": "finished"}]
        return samples

    def test_lanes_metrics(self):
        r = automatic_review.RunReview(self.run_elements_by_lane)
        lanes_metrics = r.lanes_metrics()
        assert lanes_metrics == {1: {'pc_q30': 80.7252878664961, 'yield_in_gb': 141.3262851, 'pc_pass_filter': 75.83374057009773},
                                 2: {'pc_q30': 76.55486914128899, 'yield_in_gb': 137.8495047, 'pc_pass_filter': 73.96814803233137},
                                 3: {'pc_q30': 81.09049944814811, 'yield_in_gb': 143.0777034, 'pc_pass_filter': 76.7735275382328},
                                 4: {'pc_q30': 79.0385294909938, 'yield_in_gb': 129.0695094, 'pc_pass_filter': 69.25692332762938},
                                 5: {'pc_q30': 76.01624984229616, 'yield_in_gb': 140.5316223, 'pc_pass_filter': 75.40733544260664},
                                 6: {'pc_q30': 77.22332878889476, 'yield_in_gb': 103.546344, 'pc_pass_filter': 55.561543858044104},
                                 7: {'pc_q30': 75.61055062785077, 'yield_in_gb': 143.2058529, 'pc_pass_filter': 76.84229079717159},
                                 8: {'pc_q30': 78.647185843063, 'yield_in_gb': 146.8603668, 'pc_pass_filter': 78.80325268622371}}

    def test_run_review(self):
        r = automatic_review.RunReview(self.run_elements_by_lane)
        run_review = r.run_review()
        assert run_review == {1: ['pass'],
                              2: ['pass'],
                              3: ['pass'],
                              4: ['pass'],
                              5: ['pass'],
                              6: ['fail', ['yield_in_gb']],
                              7: ['pass'],
                              8: ['pass']}

    def test_sample_metrics(self):
        s_failed = automatic_review.SampleReview(self.failed_samples)
        failed_sample_metrics = s_failed.sample_metrics()
        assert failed_sample_metrics == ('LP1251551__C_04',
                                  114.3112086,
                                  89.332777716,
                                  97.18554753953777,
                                  30.156,
                                  16.20116410973783,
                                  {'matching_snps': 28, 'no_call_chip': 1, 'mismatching_snps': 3, 'no_call_seq': 0},
                                  'male',
                                  'male')

        s_passed = automatic_review.SampleReview(self.passed_samples)
        passed_sample_metrics = s_passed.sample_metrics()
        assert passed_sample_metrics == ('LP1251553__B_06',
                                         161.185163253,
                                         127.513889915,
                                         98.22446527981002,
                                         26.716,
                                         30.956808538365465,
                                         {'no_call_seq': 0, 'no_call_chip': 2, 'mismatching_snps': 2, 'matching_snps': 28},
                                         'male',
                                         'male')

    def test_sample_review(self):
        s_failed = automatic_review.SampleReview(self.failed_samples)
        sample_review = s_failed.sample_review()
        sample_id = (list(sample_review.keys())[0])
        sample_review = (list(sample_review.values())[0])
        review_status = (sample_review[0])
        failed_fields = (sample_review[1])

        assert sample_id == 'LP1251551__C_04'
        assert review_status == 'fail'
        assert set(failed_fields) == {'clean_yield_q30', 'clean_yield_in_gb'}


