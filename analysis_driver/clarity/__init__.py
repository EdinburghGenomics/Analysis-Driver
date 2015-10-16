
from genologics.lims import Lims
from analysis_driver.config import default as cfg
from analysis_driver.app_logging import get_logger

app_logger = get_logger('clarity')

def _get_lims_connection():
    return Lims(cfg.get('clarity'))

def get_valid_lanes_from_HiseqX(flowcell_name):
    """
    Query the LIMS and return a list of valid lane for a given flowcell
    :param flowcell_name: the flowcell id such as HCH25CCXX
    :return: list of valid lane number
    """
    lims = _get_lims_connection()
    containers = lims.get_containers(type='Patterned Flowcell', name=flowcell_name)
    assert len(containers) ==1, "%s Flowcell found for name %s"%(len(containers), flowcell_name)

    flowcell = containers[0]
    valid_lanes = []
    for placement_key in flowcell.placements:
        lane = int(placement_key.split(':')[0])
        artifact = flowcell.placements.get(placement_key)
        if not artifact.udf.get('Lane Failed?', False):
            valid_lanes.append(lane)
    return valid_lanes

