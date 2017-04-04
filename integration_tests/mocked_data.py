from unittest.mock import Mock

class NamedMock(Mock):

    @property
    def name(self):
        return self.real_name


class MockedSamples(NamedMock):
    project = Mock()
    project.name = '10015AT'

mocked_lane_artifact_pool = NamedMock(
    real_name='artpool',
    input_artifact_list=Mock(return_value=[
        NamedMock(id='LP6002014-DTP_A01', reagent_labels=['D701-D502 (ATTACTCG-ATAGAGGC)'],
                  samples=[MockedSamples(real_name='10015AT0001')]),
        NamedMock(id='LP6002014-DTP_A02', reagent_labels=['D702-D502 (TCCGGAGA-ATAGAGGC)'],
                  samples=[MockedSamples(real_name='10015AT0002')]),
        NamedMock(id='LP6002014-DTP_A03', reagent_labels=['D703-D502 (CGCTCATT-ATAGAGGC)'],
                  samples=[MockedSamples(real_name='10015AT0003')]),
        NamedMock(id='LP6002014-DTP_A04', reagent_labels=['D704-D502 (GAGATTCC-ATAGAGGC)'],
                  samples=[MockedSamples(real_name='10015AT0004')]),
        NamedMock(id='LP6002014-DTP_A05', reagent_labels=['D705-D502 (ATTCAGAA-ATAGAGGC)'],
                  samples=[MockedSamples(real_name='10015AT0005')]),
        NamedMock(id='LP6002014-DTP_A06', reagent_labels=['D706-D502 (GAATTCGT-ATAGAGGC)'],
                  samples=[MockedSamples(real_name='10015AT0006')]),
        NamedMock(id='LP6002014-DTP_A07', reagent_labels=['D707-D502 (CTGAAGCT-ATAGAGGC)'],
                  samples=[MockedSamples(real_name='10015AT0007')]),
        NamedMock(id='LP6002014-DTP_A08', reagent_labels=['D708-D502 (TAATGCGC-ATAGAGGC)'],
                  samples=[MockedSamples(real_name='10015AT0008')])
    ]),
    parent_process=Mock(type=NamedMock(real_name='Create PDP Pool')),
reagent_labels=[
    'D701-D502 (ATTACTCG-ATAGAGGC)',
    'D702-D502 (TCCGGAGA-ATAGAGGC)',
    'D703-D502 (CGCTCATT-ATAGAGGC)',
    'D704-D502 (GAGATTCC-ATAGAGGC)',
    'D705-D502 (ATTCAGAA-ATAGAGGC)',
    'D706-D502 (GAATTCGT-ATAGAGGC)',
    'D707-D502 (CTGAAGCT-ATAGAGGC)',
    'D708-D502 (TAATGCGC-ATAGAGGC)'
], samples=[
        '10015AT0001',
        '10015AT0002',
        '10015AT0003',
        '10015AT0004',
        '10015AT0005',
        '10015AT0006',
        '10015AT0007',
        '10015AT0008'
    ]
)


mocked_flowcell_pooling = Mock(placements={
    '1:1': mocked_lane_artifact_pool,
    '2:1': mocked_lane_artifact_pool,
    '3:1': mocked_lane_artifact_pool,
    '4:1': mocked_lane_artifact_pool,
    '5:1': mocked_lane_artifact_pool,
    '6:1': mocked_lane_artifact_pool,
    '7:1': mocked_lane_artifact_pool,
    '8:1': mocked_lane_artifact_pool
})


class MockedRunProcess(Mock):

    def parent_processes(self):
        return [self]

    def output_containers(self):
        return [self.container]


mocked_clarity_run = MockedRunProcess(udf={'Run Status': 'RunCompleted'}, container=mocked_flowcell_pooling)
