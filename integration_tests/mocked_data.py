import os
from contextlib import contextmanager
from unittest.mock import Mock, patch
from analysis_driver.config import cfg
from analysis_driver.tool_versioning import toolset


class NamedMock(Mock):  # don't import the tests - that patches `toolset`!
    @property
    def name(self):
        return self.real_name


class MockedSamples(NamedMock):
    project = NamedMock(real_name='10015AT')


mocked_lane_artifact_pool = NamedMock(
    real_name='artpool',
    input_artifact_list=Mock(
        return_value=[
            NamedMock(reagent_labels=['D701-D502 (ATTACTCG-ATAGAGGC)'],
                      samples=[MockedSamples(real_name='10015AT0001', id='LP6002014-DTP_A01')]),
            NamedMock(reagent_labels=['D702-D502 (TCCGGAGA-ATAGAGGC)'],
                      samples=[MockedSamples(real_name='10015AT0002', id='LP6002014-DTP_A02')]),
            NamedMock(reagent_labels=['D703-D502 (CGCTCATT-ATAGAGGC)'],
                      samples=[MockedSamples(real_name='10015AT0003', id='LP6002014-DTP_A03')]),
            NamedMock(reagent_labels=['D704-D502 (GAGATTCC-ATAGAGGC)'],
                      samples=[MockedSamples(real_name='10015AT0004', id='LP6002014-DTP_A04')]),
            NamedMock(reagent_labels=['D705-D502 (ATTCAGAA-ATAGAGGC)'],
                      samples=[MockedSamples(real_name='10015AT0006', id='LP6002014-DTP_A05')]),
            NamedMock(reagent_labels=['D706-D502 (GAATTCGT-ATAGAGGC)'],
                      samples=[MockedSamples(real_name='10015AT0007', id='LP6002014-DTP_A06')]),
            NamedMock(reagent_labels=['D707-D502 (CTGAAGCT-ATAGAGGC)'],
                      samples=[MockedSamples(real_name='10015AT0008', id='LP6002014-DTP_A07')]),
            NamedMock(reagent_labels=['D708-D502 (TAATGCGC-ATAGAGGC)'],
                      samples=[MockedSamples(real_name='10015AT0009', id='LP6002014-DTP_A08')])
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


def _fake_get_list_of_samples(sample_names):
    samples = []
    for n in sample_names:
        m = Mock(udf={'Yield for Quoted Coverage (Gb)': 0.9})
        m.name = n
        samples.append(m)
    return samples


def _fake_get_user_sample_id(sample_name, lenient=False):
    return 'uid_' + sample_name


def _fake_get_plate_id_and_well(sample_name):
    return [sample_name + '_plate', 1337]


def _fake_welldups_cmd(self):
    output_file = os.path.join(self.output_directory, self.dataset.name + '.wellduplicate')
    return '{welldups} -f {coords} -r {run_dir} -t 1101 -s hiseq_x > {outfile} 2> {outfile}.err'.format(
        welldups=toolset['well_duplicates'],
        coords=cfg['well_duplicates']['coord_file'], run_dir=self.run_directory, outfile=output_file
    )


@contextmanager
def patch_pipeline(species='Homo sapiens', analysis_type='Variant Calling gatk'):
    patches = []

    def _patch(ppath, **kwargs):
        _p = patch('analysis_driver.' + ppath, **kwargs)
        _p.start()
        patches.append(_p)

    def _fake_get_sample(sample_name):
        return Mock(
            name=sample_name,
            udf={
                'Coverage': 1337,
                'Analysis Type': analysis_type,
                'Yield for Quoted Coverage (Gb)': 15,
                'Required Yield (Gb)': 30,
                'Coverage (X)': 15,
            }
        )

    _patch('client.load_config')
    _patch('dataset.clarity.get_species_from_sample', new=species)
    _patch('dataset.clarity.get_user_sample_name', new=_fake_get_user_sample_id)
    _patch('dataset.clarity.get_run', return_value=mocked_clarity_run)
    _patch('dataset.clarity.get_expected_yield_for_sample', return_value=0.9)
    _patch('dataset.clarity.get_sample', new=_fake_get_sample)
    _patch('dataset.LimsNotification')
    _patch('dataset_scanner.get_list_of_samples', new=_fake_get_list_of_samples)
    _patch('pipelines.common.clarity.find_project_name_from_sample', return_value='10015AT')
    _patch('pipelines.common.clarity.get_sample', return_value=Mock(udf={}))
    _patch('pipelines.common.clarity.get_user_sample_name', new=_fake_get_user_sample_id)
    _patch('report_generation.sample_crawler.clarity.get_species_from_sample', return_value=species)
    _patch('report_generation.sample_crawler.clarity.get_sample', new=_fake_get_sample)
    _patch('report_generation.sample_crawler.clarity.get_expected_yield_for_sample', return_value=0.9)
    _patch('report_generation.sample_crawler.clarity.get_user_sample_name', new=_fake_get_user_sample_id)
    _patch('report_generation.sample_crawler.clarity.get_plate_id_and_well', new=_fake_get_plate_id_and_well)
    _patch('report_generation.sample_crawler.clarity.get_sample_gender')
    _patch('quality_control.genotype_validation.clarity.find_project_name_from_sample', return_value='10015AT')
    _patch('quality_control.genotype_validation.clarity.get_samples_arrived_with', return_value=set())
    _patch('quality_control.genotype_validation.clarity.get_samples_genotyped_with', return_value=set())
    _patch('quality_control.genotype_validation.clarity.get_samples_sequenced_with', return_value=set())
    _patch('quality_control.genotype_validation.clarity.get_sample_names_from_project', return_value=set())
    _patch('quality_control.genotype_validation.clarity.get_sample_genotype', return_value=set())
    _patch('quality_control.well_duplicates.WellDuplicates._welldups_cmd', new=_fake_welldups_cmd)
    _patch(
        'pipelines.demultiplexing.BadTileCycleDetector',
        return_value=Mock(
            detect_bad_cycles=Mock(return_value={5: [309, 310]}),
            detect_bad_tiles=Mock(return_value={})
        )
    )

    yield

    for p in patches:
        p.stop()
