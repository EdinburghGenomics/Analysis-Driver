import os
from unittest.mock import Mock
from analysis_driver.config import cfg
from analysis_driver.tool_versioning import toolset


class NamedMock(Mock):  # don't import the tests - that patches `toolset`!
    @property
    def name(self):
        return self.real_name


mocked_clarity_project = NamedMock(real_name='10015AT', udf={'Number of Quoted Samples': 2})


class MockedSample(NamedMock):
    project = mocked_clarity_project


mocked_lane_artifact_pool = NamedMock(
    real_name='artpool',
    input_artifact_list=Mock(
        return_value=[
            NamedMock(reagent_labels=['D701-D502 (ATTACTCG-ATAGAGGC)'],
                      samples=[MockedSample(real_name='10015AT0001', id='LP6002014-DTP_A01')]),
            NamedMock(reagent_labels=['D702-D502 (TCCGGAGA-ATAGAGGC)'],
                      samples=[MockedSample(real_name='10015AT0002', id='LP6002014-DTP_A02')]),
            NamedMock(reagent_labels=['D703-D502 (CGCTCATT-ATAGAGGC)'],
                      samples=[MockedSample(real_name='10015AT0003', id='LP6002014-DTP_A03')]),
            NamedMock(reagent_labels=['D704-D502 (GAGATTCC-ATAGAGGC)'],
                      samples=[MockedSample(real_name='10015AT0004', id='LP6002014-DTP_A04')]),
            NamedMock(reagent_labels=['D705-D502 (ATTCAGAA-ATAGAGGC)'],
                      samples=[MockedSample(real_name='10015AT0006', id='LP6002014-DTP_A05')]),
            NamedMock(reagent_labels=['D706-D502 (GAATTCGT-ATAGAGGC)'],
                      samples=[MockedSample(real_name='10015AT0007', id='LP6002014-DTP_A06')]),
            NamedMock(reagent_labels=['D707-D502 (CTGAAGCT-ATAGAGGC)'],
                      samples=[MockedSample(real_name='10015AT0008', id='LP6002014-DTP_A07')]),
            NamedMock(reagent_labels=['D708-D502 (TAATGCGC-ATAGAGGC)'],
                      samples=[MockedSample(real_name='10015AT0009', id='LP6002014-DTP_A08')])
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


class MockedRunProcess(Mock):
    def parent_processes(self):
        return [self]

    def output_containers(self):
        return [self.container]


mocked_clarity_run = MockedRunProcess(udf={'Run Status': 'RunCompleted'}, container=mocked_flowcell_pooling)


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
