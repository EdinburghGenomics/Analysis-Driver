import os
from unittest.mock import Mock
from analysis_driver.config import cfg
from analysis_driver.tool_versioning import toolset
from tests.test_analysisdriver import NamedMock


mocked_clarity_project = NamedMock(real_name='10015AT', udf={'Number of Quoted Samples': 2})


class MockedSample(NamedMock):
    project = mocked_clarity_project


mocked_lane_artifact_pool = NamedMock(
    real_name='artpool',
    input_artifact_list=Mock(
        return_value=[
            Mock(reagent_labels=['D701-D502 (ATTACTCG-ATAGAGGC)'], samples=[MockedSample(real_name='10015AT0001', id='LP6002014-DTP_A01')]),
            Mock(reagent_labels=['D702-D502 (TCCGGAGA-ATAGAGGC)'], samples=[MockedSample(real_name='10015AT0002', id='LP6002014-DTP_A02')]),
            Mock(reagent_labels=['D703-D502 (CGCTCATT-ATAGAGGC)'], samples=[MockedSample(real_name='10015AT0003', id='LP6002014-DTP_A03')]),
            Mock(reagent_labels=['001A IDT-ILMN TruSeq DNA-RNA UD 96 Indexes  Plate_UDI0001 (GAGATTCC-ATAGAGGC)'], samples=[MockedSample(real_name='10015AT0004', id='LP6002014-DTP_A04')]),
            Mock(reagent_labels=['D705-D502 (ATTCAGAA-ATAGAGGC)'], samples=[MockedSample(real_name='10015AT0006', id='LP6002014-DTP_A05')]),
            Mock(reagent_labels=['D706-D502 (GAATTCGT-ATAGAGGC)'], samples=[MockedSample(real_name='10015AT0007', id='LP6002014-DTP_A06')]),
            Mock(reagent_labels=['D707-D502 (CTGAAGCT-ATAGAGGC)'], samples=[MockedSample(real_name='10015AT0008', id='LP6002014-DTP_A07')]),
            Mock(reagent_labels=['D708-D502 (TAATGCGC-ATAGAGGC)'], samples=[MockedSample(real_name='10015AT0009', id='LP6002014-DTP_A08')])
        ]
    ),
    parent_process=Mock(type=NamedMock(real_name='Create PDP Pool')),
    reagent_labels=['D701-D502 (ATTACTCG-ATAGAGGC)', 'D702-D502 (TCCGGAGA-ATAGAGGC)',
                    'D703-D502 (CGCTCATT-ATAGAGGC)', 'D704-D502 (GAGATTCC-ATAGAGGC)',
                    'D705-D502 (ATTCAGAA-ATAGAGGC)', 'D706-D502 (GAATTCGT-ATAGAGGC)',
                    'D707-D502 (CTGAAGCT-ATAGAGGC)', 'D708-D502 (TAATGCGC-ATAGAGGC)'],
    samples=['10015AT0001', '10015AT0002', '10015AT0003', '10015AT0004', '10015AT0005', '10015AT0006',
             '10015AT0007', '10015AT0008']
)


mocked_flowcell_pooling = Mock(
    placements={'1:1': mocked_lane_artifact_pool, '2:1': mocked_lane_artifact_pool,
                '3:1': mocked_lane_artifact_pool, '4:1': mocked_lane_artifact_pool,
                '5:1': mocked_lane_artifact_pool, '6:1': mocked_lane_artifact_pool,
                '7:1': mocked_lane_artifact_pool, '8:1': mocked_lane_artifact_pool}
)


def fake_non_pooling_artifact(i, rapid='no'):
    return Mock(
        reagent_labels=['A001-B002 (ATGCATGC-CTGACTGA)'],
        samples=[
            MockedSample(
                real_name='non_pooling_sample_' + i,
                id='a_library',
                udf={'Rapid Analysis': rapid, 'User Sample Name': 'uid_non_pooling_sample_' + i}
            )
        ]
    )

mocked_flowcell_non_pooling = Mock(
    placements={
        '1:1': fake_non_pooling_artifact('1', 'No'),
        '2:1': fake_non_pooling_artifact('2', 'Yes'),
        '3:1': fake_non_pooling_artifact('3', 'No'),
        '4:1': fake_non_pooling_artifact('4', 'No'),
        '5:1': fake_non_pooling_artifact('5', 'No'),
        '6:1': fake_non_pooling_artifact('6', 'No'),
        '7:1': fake_non_pooling_artifact('7', 'No'),
        '8:1': fake_non_pooling_artifact('8', 'Yes')
    }
)


class MockedRunProcess(Mock):
    def parent_processes(self):
        return [self]

    def output_containers(self):
        return [self.container]


mocked_pooling_run = MockedRunProcess(udf={'Run Status': 'RunCompleted'}, container=mocked_flowcell_pooling)
mocked_rapid_run = MockedRunProcess(udf={'Run Status': 'RunCompleted'}, container=mocked_flowcell_non_pooling)


def fake_get_user_sample_id(sample_name, lenient=False):
    return 'uid_' + sample_name


def fake_get_plate_id_and_well(sample_name):
    return [sample_name + '_plate', 1337]


def fake_welldups_cmd(self):
    output_file = os.path.join(self.output_directory, self.dataset.name + '.wellduplicate')
    return '{welldups} -f {coords} -r {run_dir} -t 1101 -s hiseq_x > {outfile} 2> {outfile}.err'.format(
        welldups=toolset['well_duplicates'],
        coords=cfg['well_duplicates']['coord_file'], run_dir=self.run_directory, outfile=output_file
    )

sample_1 = {
    'project_id': "project1",
    'sample_id': "sample1",
    'barcode': 'D701-D502 (ATTACTCG-ATAGAGGC)',
    'artifact_id': "2-1111111"
}

sample_2 = {
    'project_id': "project1",
    'sample_id': "sample2",
    'barcode': 'D702-D502 (TCCGGAGA-ATAGAGGC)',
    'artifact_id': "2-1111112"
}

sample_3 = {
    'project_id': "project1",
    'sample_id': "sample3",
    'barcode': 'D703-D502 (CGCTCATT-ATAGAGGC)',
    'artifact_id': "2-1111113"
}

sample_4 = {
    'project_id': "project1",
    'sample_id': "sample4",
    'barcode': "002H IDT-ILMN TruSeq DNA-RNA UD 96 Indexes Plate_UDI0016 (AATCCGGA-CTACAGTT)",
    'artifact_id': "2-1111114"
}

sample_5 = {
    'project_id': "project1",
    'sample_id': "sample5",
    'barcode': 'D707-D502 (CTGAAGCT-ATAGAGGC)',
    'artifact_id': "2-1111115"
}

sample_6 = {
    'project_id': "project1",
    'sample_id': "sample6",
    'barcode': 'D707-D502 (CTGAAGCT-ATAGAGGC)',
    'artifact_id': "2-1111116"
}

sample_7 = {
    'project_id': "project1",
    'sample_id': "sample7",
    'barcode': 'D708-D502 (TAATGCGC-ATAGAGGC)',
    'artifact_id': "2-1111117"
}

sample_8 = {
    'project_id': "project1",
    'sample_id': "sample8",
    'barcode': 'D708-D502 (TAATGCGC-ATAGAGGC)',
    'artifact_id': "2-1111118"
}

pool_of_samples = [sample_1, sample_2, sample_3, sample_4]

mocked_run_status = {
    'run_status': "RunCompleted",
    'lanes': [
        {'lane': 1, 'samples': [sample_1]},
        {'lane': 2, 'samples': [sample_2]},
        {'lane': 3, 'samples': [sample_3]},
        {'lane': 4, 'samples': [sample_4]},
        {'lane': 5, 'samples': [sample_5]},
        {'lane': 6, 'samples': [sample_6]},
        {'lane': 7, 'samples': [sample_7]},
        {'lane': 8, 'samples': [sample_8]},
    ],
    'instrument_id': "HSX-E00375",
    'nb_reads': "3",
    'nb_cycles': "310"
}

mocked_run_status_pools = {
    'run_status': "RunCompleted",
    'lanes': [
        {'lane': 1, 'samples': pool_of_samples},
        {'lane': 2, 'samples': pool_of_samples},
        {'lane': 3, 'samples': pool_of_samples},
        {'lane': 4, 'samples': pool_of_samples},
        {'lane': 5, 'samples': pool_of_samples},
        {'lane': 6, 'samples': pool_of_samples},
        {'lane': 7, 'samples': pool_of_samples},
        {'lane': 8, 'samples': pool_of_samples},
    ],
    'instrument_id': "HSX-E00375",
    'nb_reads': "3",
    'nb_cycles': "310"
}

