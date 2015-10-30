import re
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


def sanitize_user_id(user_id):
    if isinstance(str,user_id):
        return re.sub("[^\w_\-.]","_",user_id)
    else:
        return None


def get_user_sample_name(sample_name):
    """
    Query the LIMS and return the name the user gave to the sample
    :param sample_name: the sample name from the Samplesheet.csv
    :return: the user's sample name or None
    """
    lims = _get_lims_connection()
    samples = lims.get_samples(name=sample_name)
    #FIXME: Remove the hack when we're sure our sample id don't have colon
    if len(samples) == 0:
        sample_name_sub = re.sub("_(\d{2})",":\g<1>",sample_name)
        samples = lims.get_samples(name=sample_name_sub)
    if len(samples) != 1:
        app_logger.warning('%s Sample(s) found for name %s' % (len(samples), sample_name))
        return None

    return sanitize_user_id(samples[0].udf.get('User Sample Name'))


def run_tests():
    assert get_valid_lanes('HCH25CCXX') == [1, 2, 3, 4, 5, 6, 7]
    assert get_valid_lanes('HCH25CCX') is None

    assert get_user_sample_name('10094AT0001') == '1118-RP'
    assert get_user_sample_name('NA12877_25SEPT15 2/5') is None

if __name__ == '__main__':
    # will only work with a valid connection to the production server
    run_tests()
