import sqlite3
from unittest.mock import patch
from tests.test_external_data import FakeRestResponse
from analysis_driver.external_data import ncbi


def reset_cache():
    ncbi.data_cache = sqlite3.connect(':memory:')
    ncbi.cursor = ncbi.data_cache.cursor()
    ncbi.cursor.execute('CREATE TABLE IF NOT EXISTS species (taxid text UNIQUE, scientific_name text UNIQUE, common_name text)')
    ncbi.cursor.execute('CREATE TABLE IF NOT EXISTS aliases (query_name text UNIQUE, taxid text REFERENCES species(taxid))')


def test_fetch_from_eutils():
    ncbi_search_data = {'esearchresult': {'idlist': ['1337']}}
    ncbi_fetch_data = '''
        <ScientificName>Genus species</ScientificName>
        <OtherNames><CommonName>a common name</CommonName></OtherNames>
        <Rank>species</Rank>
    '''

    patched_get = patch(
        'analysis_driver.external_data.ncbi.requests.get',
        side_effect=(
            FakeRestResponse(content=ncbi_search_data),
            FakeRestResponse(content=ncbi_fetch_data),
            FakeRestResponse(content=ncbi_fetch_data),
            FakeRestResponse(content=ncbi_fetch_data)
        )
    )

    with patched_get as mocked_get:
        obs = ncbi._fetch_from_eutils('a_species')
        assert obs == ('1337', 'Genus species', 'a common name')
        mocked_get.assert_any_call(
            'http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi',
            params={'db': 'Taxonomy', 'term': 'a_species', 'retmode': 'JSON'}
        )
        mocked_get.assert_any_call(
            'http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi',
            params={'db': 'Taxonomy', 'id': '1337'}
        )


def test_cache():
    assert ncbi._fetch_from_cache('a species') is None
    ncbi._cache_species('a species', '1337', 'Scientific name', 'a species')
    assert ncbi._fetch_from_cache('a species') == ('a species', '1337', 'Scientific name', 'a species')
    reset_cache()


def test_get_species_name():
    fetch = 'analysis_driver.external_data.ncbi._fetch_from_eutils'
    assert ncbi._fetch_from_cache('a species') is None
    with patch(fetch, return_value=(None, None, None)):
        assert ncbi.get_species_name('a species') is None
        assert ncbi._fetch_from_cache('a species') is None
    reset_cache()
    with patch(fetch, return_value=('1337', 'Scientific name', 'a species')):
        assert ncbi.get_species_name('a species') == 'Scientific name'
        assert ncbi._fetch_from_cache('a species') == ('a species', '1337', 'Scientific name', 'a species')
    reset_cache()
