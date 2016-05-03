from unittest.mock import patch, Mock
import os
from analysis_driver import clarity
from tests.test_analysisdriver import TestAnalysisDriver
from tests.test_rest_communication import FakeRestResponse

clarity._lims = Mock()
helper = TestAnalysisDriver()


def patched(path, **kwargs):
    return patch('analysis_driver.clarity.' + path, **kwargs)


def patched_lims(method, return_value=None, side_effect=None):
    return patched('_lims.' + method, return_value=return_value, side_effect=side_effect)


def patched_clarity(function, return_value=None, side_effect=None):
    return patched(function, return_value=return_value, side_effect=side_effect)


class FakeEntity(Mock):
    def __init__(self, name, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.name = name


class FakeContainer:
    @staticmethod
    def get_placements():
        return {
            'this': Mock(samples=[FakeEntity('a name')]),
            'that': Mock(samples=[FakeEntity('another name')])
        }


class FakeProcess:
    @staticmethod
    def all_inputs():
        return [Mock(samples=[FakeEntity('this'), FakeEntity('that')])]

    @staticmethod
    def input_per_sample(sample_name):
        return [Mock(id=sample_name)]

    @staticmethod
    def outputs_per_input(artifact_id, **kwargs):
        return [Mock(container=artifact_id)]


fake_samples = [
    Mock(project=FakeEntity('this'), udf={'Species': 'a_species'}),
    Mock(project=FakeEntity('this'), udf={'Species': 'a_species'}),
    Mock(project=FakeEntity('that'), udf={'Species': 'a_species'})
]


def test_get_valid_lanes():
    fake_flowcell = Mock(
        placements={
            '1:this': Mock(udf={'Lane Failed?': False}),
            '2:that': Mock(udf={'Lane Failed?': False}),
            '3:other': Mock(udf={'Lane Failed?': True})
        }
    )
    with patched_lims('get_containers', [fake_flowcell]) as mocked_lims:
        valid_lanes = clarity.get_valid_lanes('a_flowcell_name')
        mocked_lims.assert_called_with(type='Patterned Flowcell', name='a_flowcell_name')
        assert valid_lanes == [1, 2]


def test_find_project_from_sample():
    with patched_clarity('get_samples', fake_samples) as mocked_get_samples:
        project_name = clarity.find_project_name_from_sample('a_sample')
        mocked_get_samples.assert_called_with('a_sample')
        assert project_name is None

    with patched_clarity('get_samples', fake_samples[0:1]):
        assert clarity.find_project_name_from_sample('a_sample') == 'this'


@patched_lims('get_artifacts', [Mock(parent_process=Mock(udf={}, input_per_sample=lambda sample_name: [Mock(position='1:this', udf={})]))])
@patched_clarity('get_sample', FakeEntity('a_sample'))
def test_find_run_elements_from_sample(mocked_get_sample, mocked_get_artifacts):
    assert list(clarity.find_run_elements_from_sample('a_sample')) == [(None, '1')]
    mocked_get_sample.assert_called_with('a_sample')
    mocked_get_artifacts.assert_called_with(sample_name='a_sample', process_type='AUTOMATED - Sequence')


def test_get_species_information_from_ncbi():
    ncbi_search_data = {'esearchresult': {'idlist': ['1337']}}
    ncbi_fetch_data = '''
        <ScientificName>Genus species</ScientificName>
        <OtherNames><CommonName>a common name</CommonName></OtherNames>
        <Rank>species</Rank>
    '''

    patched_get = patch(
        'analysis_driver.clarity.requests.get',
        side_effect=(
            FakeRestResponse(content=ncbi_search_data),
            FakeRestResponse(content=ncbi_fetch_data),
            FakeRestResponse(content=ncbi_fetch_data),
            FakeRestResponse(content=ncbi_fetch_data)
        )
    )

    with patched_get as mocked_get:
        obs = clarity.get_species_information_from_ncbi('a_species')
        assert obs == ('1337', 'Genus species', 'a common name')
        mocked_get.assert_any_call(
            'http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi',
            params={'db': 'Taxonomy', 'term': 'a_species', 'retmode': 'JSON'}
        )
        mocked_get.assert_any_call(
            'http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi',
            params={'db': 'Taxonomy', 'id': '1337'}
        )


@patched_clarity('get_species_information_from_ncbi', ('1337', 'Genus species', 'common name'))
@patched_clarity('get_samples', fake_samples)
def test_get_species_from_sample(mocked_get_samples, mocked_ncbi):
    assert clarity.get_species_from_sample('a_sample_name') == 'Genus species'
    mocked_get_samples.assert_called_with('a_sample_name')
    mocked_ncbi.assert_called_with('a_species')


def test_sanitize_user_id():
    assert clarity.sanitize_user_id('this?that$other another:more') == 'this_that_other_another_more'


def test_get_list_of_samples():
    exp_lims_sample_ids = ['this', 'that:01', 'other _L:01']
    calling_sample_ids = ['this', 'that_01', 'other__L_01']
    fake_list_samples = [[FakeEntity(n)] for n in exp_lims_sample_ids]
    pbatch = patched_lims('get_batch')
    psamples = patched_lims('get_samples', side_effect=fake_list_samples)

    with pbatch, psamples as mocked_get_samples:
        samples = clarity.get_list_of_samples(calling_sample_ids)
        assert [s.name for s in samples] == exp_lims_sample_ids
        mocked_get_samples.assert_any_call(name=['this', 'that_01', 'other__L_01'])
        mocked_get_samples.assert_any_call(name=['other__L:01', 'that:01'])
        mocked_get_samples.assert_any_call(name=['other _L:01'])


@patched_lims('get_samples', side_effect=[[], [], [None]])
def test_get_samples(mocked_lims):
    assert clarity.get_samples('a_sample_name__L_01') == [None]
    mocked_lims.assert_any_call(name='a_sample_name__L_01')
    mocked_lims.assert_any_call(name='a_sample_name__L:01')
    mocked_lims.assert_any_call(name='a_sample_name _L:01')


@patched_clarity('get_samples', return_value=['a sample'])
def test_get_sample(mocked_lims):
    assert clarity.get_sample('a_sample_id') == 'a sample'
    mocked_lims.assert_called_with('a_sample_id')


@patched_clarity('get_sample', return_value=Mock(udf={'User Sample Name': 'a:user:sample:id'}))
def test_get_user_sample_name(mocked_lims):
    assert clarity.get_user_sample_name('a_sample_id') == 'a_user_sample_id'
    mocked_lims.assert_called_with('a_sample_id')


@patched_clarity('get_sample', return_value=Mock(udf={'Gender': 'unknown'}))
def test_get_sample_gender(mocked_lims):
    assert clarity.get_sample_gender('a_sample_id') == 'unknown'
    mocked_lims.assert_called_with('a_sample_id')


@patched_lims('get_file_contents', 'some test content')
@patched_clarity('get_sample', Mock(udf={'Genotyping results file id': 1337}))
def test_get_genotype_information_from_lims(mocked_get_sample, mocked_file_contents):
    genotype_vcf = os.path.join(helper.assets_path, 'a_genotype.vcf')
    assert clarity.get_sample_genotype('a_sample_name', genotype_vcf) == genotype_vcf
    mocked_get_sample.assert_called_with('a_sample_name')
    mocked_file_contents.assert_called_with(id=1337)
    assert open(genotype_vcf).read() == 'some test content'
    os.remove(genotype_vcf)


@patched_clarity('get_sample', Mock(udf={'Yield for Quoted Coverage (Gb)': 3}))
def test_get_expected_yield_for_sample(mocked_get_sample):
    assert clarity.get_expected_yield_for_sample('a_sample_id') == 3000000000
    mocked_get_sample.assert_called_with('a_sample_id')


@patched_lims('get_processes', ['a_run'])
def test_get_run(mocked_lims):
    assert clarity.get_run('a_run_id') == 'a_run'
    mocked_lims.assert_called_with(type='AUTOMATED - Sequence', udf={'RunID': 'a_run_id'})


@patched_lims('route_artifacts')
@patched_clarity('get_sample', side_effect=[Mock(artifact='this'), Mock(artifact='that')])
@patched_lims('get_uri', 'a_workflow_uri')
def test_route_samples_to_delivery_workflow(mocked_get_uri, mocked_get_sample, mocked_route):
    clarity.route_samples_to_delivery_workflow(['a_sample_id', 'another_sample_id'])
    mocked_get_uri.assert_called_with('configuration', 'workflows', '401')
    mocked_get_sample.assert_any_call('a_sample_id')
    mocked_get_sample.assert_any_call('another_sample_id')
    mocked_route.assert_called_with(['this', 'that'], workflow_uri='a_workflow_uri')


@patched_clarity('get_samples', [Mock(artifact=Mock(location=(FakeEntity('a_plate'), 'a_well')))])
def test_get_plate_id_and_well_from_lims(mocked_lims):
    assert clarity.get_plate_id_and_well('a_sample_id') == ('a_plate', 'a_well')
    mocked_lims.assert_called_with('a_sample_id')


@patched_lims('get_containers', [FakeContainer])
def test_get_sample_names_from_plate_from_lims(mocked_lims):
    obs = clarity.get_sample_names_from_plate('a_plate_id')
    assert sorted(obs) == ['a_name', 'another_name']
    mocked_lims.assert_called_with(type='96 well plate', name='a_plate_id')


@patched_lims('get_samples', [FakeEntity('this'), FakeEntity('that')])
def test_get_sample_names_from_project_from_lims(mocked_lims):
    assert clarity.get_sample_names_from_project('a_project') == ['this', 'that']
    mocked_lims.assert_called_with(projectname='a_project')


@patched_lims('get_processes', [FakeProcess])
@patched_lims('get_artifacts', [Mock(id='this'), Mock(id='that')])
@patched_clarity('get_sample', FakeEntity('a_sample_name'))
def test_get_output_containers_from_sample_and_step_name(mocked_get_sample, mocked_get_arts, mocked_get_prcs):
    obs = clarity.get_output_containers_from_sample_and_step_name('a_sample_id', 'a_step_name')
    assert obs == {'a_sample_name'}
    mocked_get_sample.assert_called_with('a_sample_id')
    mocked_get_arts.assert_called_with(sample_name='a_sample_name')
    mocked_get_prcs.assert_called_with(type='a_step_name', inputartifactlimsid=['this', 'that'])


@patched_clarity('get_sample_names_from_plate', ['this', 'that', 'other'])
@patched_clarity('get_sample',
                 Mock(artifact=Mock(container=FakeEntity('a_container', type=FakeEntity('96 well plate')))))
def test_get_samples_arrived_with(mocked_get_sample, mocked_names_from_plate):
    assert clarity.get_samples_arrived_with('a_sample_name') == ['this', 'that', 'other']
    mocked_get_sample.assert_called_with('a_sample_name')
    mocked_names_from_plate.assert_called_with('a_container')


@patched_clarity('get_sample_names_from_plate', ['other'])
@patched_clarity('get_output_containers_from_sample_and_step_name', [FakeEntity('this'), FakeEntity('that')])
@patched_clarity('get_sample', FakeEntity('a_sample_name'))
def test_get_samples_genotyped_with(mocked_get_sample, mocked_containers, mocked_names_from_plate):
    assert clarity.get_samples_genotyped_with('a_sample_name') == {'other'}
    mocked_get_sample.assert_called_with('a_sample_name')
    mocked_containers.assert_called_with('a_sample_name', 'Genotyping Plate Preparation EG 1.0')
    mocked_names_from_plate.assert_any_call('this')
    mocked_names_from_plate.assert_any_call('that')


@patched_clarity('get_sample_names_from_plate', ['other'])
@patched_clarity('get_output_containers_from_sample_and_step_name', [FakeEntity('this'), FakeEntity('that')])
@patched_clarity('get_sample', FakeEntity('a_sample_name'))
def test_get_samples_sequenced_with(mocked_get_sample, mocked_containers, mocked_names_from_plate):
    assert clarity.get_samples_sequenced_with('a_sample_name') == {'other'}
    mocked_get_sample.assert_called_with('a_sample_name')
    mocked_containers.assert_called_with('a_sample_name', 'Sequencing Plate Preparation EG 1.0')
    mocked_names_from_plate.assert_any_call('this')
    mocked_names_from_plate.assert_any_call('that')


@patched_lims('get_processes', [FakeProcess])
def test_get_released_samples(mocked_lims):
    assert clarity.get_released_samples() == ['that', 'this']
    mocked_lims.assert_called_with(type='Data Release EG 1.0')
