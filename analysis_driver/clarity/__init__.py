import re
import requests
from genologics.lims import Lims
from analysis_driver.config import default as cfg
from analysis_driver.app_logging import get_logger

app_logger = get_logger('Clarity')


def _get_lims_connection():
    return Lims(**cfg.get('clarity'))


def get_valid_lanes(flowcell_name):
    """
    Query the LIMS and return a list of valid lane for a given flowcell
    :param flowcell_name: the flowcell id such as HCH25CCXX
    :return: list of valid lane number
    """
    lims = _get_lims_connection()
    containers = lims.get_containers(type='Patterned Flowcell', name=flowcell_name)
    if len(containers) != 1:
        app_logger.warning('%s Flowcell(s) found for name %s' % (len(containers), flowcell_name))
        return None

    flowcell = containers[0]
    valid_lanes = []
    for placement_key in flowcell.placements:
        lane = int(placement_key.split(':')[0])
        artifact = flowcell.placements.get(placement_key)
        if not artifact.udf.get('Lane Failed?', False):
            valid_lanes.append(lane)
    valid_lanes = sorted(valid_lanes)
    app_logger.info('Valid lanes for %s: %s' % (flowcell_name, str(valid_lanes)))
    return valid_lanes


def find_project_from_sample(sample_name):
    """Query clarity to get the project name of a sample"""
    lims = _get_lims_connection()
    samples = get_lims_samples(sample_name, lims)
    if samples:
        project_names = set([s.project.name for s in samples])
        if len(project_names) != 1:
            app_logger.error('%s projects found for sample %s' % (len(project_names), sample_name))
            return None
        else:
            return project_names.pop()


def find_run_elements_from_sample(sample_name):
    lims = _get_lims_connection()
    sample = get_lims_sample(sample_name, lims)
    if sample:
        run_log_files = lims.get_artifacts(sample_name=sample.name, process_type="AUTOMATED - Sequence")
        for run_log_file in run_log_files:
            p = run_log_file.parent_process
            run_id = p.udf.get('RunID')
            lanes = p.input_per_sample(sample.name)
            for artifact in lanes:
                lane = artifact.position.split(':')[0]
                if not artifact.udf.get('Lane Failed?', False):
                    yield run_id, lane

def get_species_information_from_ncbi(species):
    """Query NCBI taxomomy database to get the taxomoy id scientific name and common name
    Documentation available at http://www.ncbi.nlm.nih.gov/books/NBK25499/"""
    esearch_url = 'http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi'
    payload = {'db': 'Taxonomy', 'term': species, 'retmode': 'JSON'}
    r = requests.get(esearch_url, params=payload)
    results = r.json()
    taxid_list = taxid = results.get('esearchresult').get('idlist')
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
    if len(all_species_names) > 1:
        app_logger.error("More than one taxon corresponding to %s" % species)
        return None, None, None
    elif len(all_species_names) == 0:
        app_logger.error("No taxon found corresponding to %s" % species)
        return None, None, None

    return all_species_names[0]

def get_species_from_sample(sample_name):
    lims = _get_lims_connection()
    samples = get_lims_samples(sample_name, lims)
    species_string = None
    if samples:
        species_strings = set([s.udf.get('Species') for s in samples])
        if len(species_strings) != 1:
            app_logger.error('%s species found for sample %s' % (len(species_strings), sample_name))
        else:
            species_string = species_strings.pop()
    if species_string:
        taxid, scientific_name, common_name = get_species_information_from_ncbi(species_string)
        if taxid:
            return scientific_name
    return None

def sanitize_user_id(user_id):
    if isinstance(user_id, str):
        return re.sub("[^\w_\-.]", "_", user_id)
    else:
        return None


def get_lims_samples(sample_name, lims):
    samples = lims.get_samples(name=sample_name)
    # FIXME: Remove the hack when we're sure our sample id don't have colon
    if len(samples) == 0:
        sample_name_sub = re.sub("_(\d{2})$", ":\g<1>", sample_name)
        samples = lims.get_samples(name=sample_name_sub)
    if len(samples) == 0:
        sample_name_sub = re.sub("__(\w)_(\d{2})", " _\g<1>:\g<2>", sample_name)
        samples = lims.get_samples(name=sample_name_sub)
    return samples


def get_lims_sample(sample_name, lims):
    samples = get_lims_samples(sample_name, lims)
    if len(samples) != 1:
        app_logger.warning('%s Sample(s) found for name %s' % (len(samples), sample_name))
        return None
    return samples[0]


def get_user_sample_name(sample_name):
    """
    Query the LIMS and return the name the user gave to the sample
    :param sample_name: the sample name from the Samplesheet.csv
    :return: the user's sample name or None
    """
    lims = _get_lims_connection()
    sample = get_lims_sample(sample_name, lims)
    if sample:
        return sanitize_user_id(sample.udf.get('User Sample Name'))


def get_sex_from_lims(sample_name):
    lims = _get_lims_connection()
    samples = get_lims_samples(sample_name, lims)
    if len(samples) == 1:
        gender = samples[0].udf.get('Gender')
        return gender


def get_genotype_information_from_lims(sample_name, output_file_name):
    lims = _get_lims_connection()
    sample = get_lims_sample(sample_name, lims)
    if sample:
        file_id = sample.udf('genotype_results')
        if file_id:
            file_content = lims.get_file_contents(id=file_id)
            with open(file_content, 'w') as open_file:
                open_file.write(file_content)
        return output_file_name
    return None



def get_expected_yield_for_sample(sample_name):
    """
    Query the LIMS and return the number of bases expected for a sample
    :param sample_name: the sample name
    :return: number of bases
    """
    lims = _get_lims_connection()
    sample = get_lims_sample(sample_name, lims)
    if sample:
        nb_gb = sample.udf.get('Yield for Quoted Coverage (Gb)')
        if nb_gb:
            return nb_gb * 1000000000



def run_tests():
    assert get_valid_lanes('HCH25CCXX') == [1, 2, 3, 4, 5, 6, 7]
    assert get_valid_lanes('HCH25CCX') is None

    assert get_user_sample_name('10094AT0001') == '1118-RP'
    assert get_user_sample_name('NA12877_25SEPT15 2/5') is None

    assert find_run_elements_from_sample('10094AT0001')
    assert get_species_from_sample('10094AT0001') == "Homo sapiens"
    print(get_expected_yield_for_sample('X0002DM003_A_06'))
    print(get_species_information_from_ncbi('mouse'))
    print(get_species_information_from_ncbi('human'))
    print(get_species_information_from_ncbi('pig'))
    print(get_species_information_from_ncbi('fruit fly'))
    print(get_species_information_from_ncbi('chicken'))
    print(get_species_information_from_ncbi('potato'))
    print(get_species_information_from_ncbi('wheat'))


if __name__ == '__main__':
    # will only work with a valid connection to the production server
    run_tests()
