import re
import sqlite3
import requests
from analysis_driver.config import default as cfg
from analysis_driver.app_logging import logging_default as log_cfg

app_logger = log_cfg.get_logger('ncbi')


data_cache = sqlite3.connect(cfg['ncbi_cache'])
cursor = data_cache.cursor()
cursor.execute('CREATE TABLE IF NOT EXISTS species (taxid text UNIQUE, scientific_name text UNIQUE, common_name text)')
cursor.execute('CREATE TABLE IF NOT EXISTS aliases (query_name text UNIQUE, taxid text REFERENCES species(taxid))')


def get_species_name(query_species):
    local_query = _fetch_from_cache(query_species)
    if local_query:
        q, taxid, scientific_name, common_name = local_query
    else:
        taxid, scientific_name, common_name = _fetch_from_eutils(query_species)
        _cache_species(query_species, taxid, scientific_name, common_name)

    return scientific_name


def _fetch_from_cache(query_species):
    cursor.execute('SELECT * FROM aliases NATURAL JOIN species WHERE query_name=?', (query_species,))
    return cursor.fetchone()


def _cache_species(query_species, taxid, scientific_name, common_name):
    cursor.execute('SELECT taxid FROM species WHERE taxid=?', (taxid,))
    if not cursor.fetchone():
        cursor.execute('INSERT INTO species VALUES (?, ?, ?)', (taxid, scientific_name, common_name))
    cursor.execute('INSERT INTO aliases VALUES (?, ?)', (query_species, taxid))
    data_cache.commit()


def _fetch_from_eutils(species):
    """
    Query NCBI taxomomy database to get the taxomoy id scientific name and common name
    Documentation available at http://www.ncbi.nlm.nih.gov/books/NBK25499/
    :param str species: A search term for esearch, e.g. Human, Mouse, etc.
    """
    esearch_url = 'http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi'
    payload = {'db': 'Taxonomy', 'term': species, 'retmode': 'JSON'}
    r = requests.get(esearch_url, params=payload)
    results = r.json()
    taxid_list = results.get('esearchresult').get('idlist')
    all_species_names = []
    for taxid in taxid_list:
        efetch_url = 'http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi'
        payload = {'db': 'Taxonomy', 'id': taxid}
        r = requests.get(efetch_url, params=payload)
        match = re.search('<Rank>(.+?)</Rank>', r.text, re.MULTILINE)
        rank = None
        if match:
            rank = match.group(1)
        if rank == 'species':
            scientific_name = common_name = None
            match = re.search('<ScientificName>(.+?)</ScientificName>', r.text, re.MULTILINE)
            if match:
                scientific_name = match.group(1)
            match = re.search('<GenbankCommonName>(.+?)</GenbankCommonName>', r.text, re.MULTILINE)
            if not match:
                match = re.search('<CommonName>(.+?)</CommonName>', r.text, re.MULTILINE)
            if match:
                common_name = match.group(1)
            all_species_names.append((taxid, scientific_name, common_name))
    nspecies = len(all_species_names)
    if nspecies != 1:
        app_logger.error('%s taxons found for %s', nspecies, species)
        return None, None, None

    return all_species_names[0]
