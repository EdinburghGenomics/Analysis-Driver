import glob
import os
import shutil
from xml.etree import ElementTree as et

from egcg_core.constants import ELEMENT_PROJECT_ID, ELEMENT_LANE, ELEMENT_SAMPLE_INTERNAL_ID


def get_run_element_id(run_id, lane_number, barcode=None):
    """Retrieving run_element_id, which differs between barcodeless and multiplexed runs"""
    if barcode is None:
        return '{dataset_name}_{lane}'.format(dataset_name=run_id, lane=lane_number)
    else:
        return '{dataset_name}_{lane}_{barcode}'.format(dataset_name=run_id, lane=lane_number, barcode=barcode)


def prepend_path_to_data_files(prepend_path, data_structure):
    """
    Iterate over a dictionary containing relative file paths and prepend a path then store in a new data structure.
    It recurses if one of the value is found to be a dict.
    This function does not modify the provided data_structure.
    """
    new_data_structure = {}
    for k in data_structure:
        if isinstance(data_structure[k], dict):
            new_data_structure[k] = prepend_path_to_data_files(prepend_path, data_structure[k])
        else:
            new_data_structure[k] = os.path.join(prepend_path, data_structure[k])
    return new_data_structure


def split_in_chunks(total_length, chunksize, zero_based=False, end_inclusive=True):
    """
    Create  N chunks where:
      -  N is the number of full chunks (nb_lines // split_lines)
      -  + 1 or 0 if there are no extra incomplete chunks
    The indices can be 0 or 1 based but will all be end inclusive
    :param total_length: The length to split in chunks
    :param chunksize: The size of the chunk
    :param zero_based: True if the start position where to start the first chunk
    :param end_inclusive: True if the end of the chunk is included False otherwise
    :return: List of tuples where first element is start and second is the end of the chunk created.
    """
    last = 0 if zero_based else 1
    chunks = []
    chunk_padding = 1 if end_inclusive else 0
    if end_inclusive and zero_based:
        end_padding = -1
    elif end_inclusive or zero_based:
        end_padding = 0
    else:
        end_padding = 1
    for chunki in range(total_length // chunksize + min(1, total_length % chunksize)):
        new = last + chunksize - chunk_padding
        chunks.append((last, min(new, total_length + end_padding)))
        last = new + chunk_padding
    return chunks


def merge_lane_directories(fastq_dir, run_elements):
    """
    At the end of a per lane bcl2fastq run the data in the lane directories must be merged back to the main fastq dir
    This function does that by moving all the fastq files to the expected location
    and the stats file are move to a central Stats folder. They are rename to avoid file collision.
    :param fastq_dir: The main fastq dir
    :param run_elements: The run element expected in this run
    """
    # find all the lane directories:
    lane_dirs = glob.glob(os.path.join(fastq_dir, 'lane_*'))

    os.makedirs(os.path.join(fastq_dir, 'Stats'), exist_ok=True)

    # find the stats files and unassigned if they exist
    for lane_dir in lane_dirs:
        # We're working with lane
        lane_num = lane_dir.split('_')[-1]
        # Move each file in the stats directory
        stats_dir = os.path.join(lane_dir, 'Stats')
        for fname in os.listdir(stats_dir):
            os.rename(
                os.path.join(stats_dir, fname),
                os.path.join(fastq_dir, 'Stats', 'lane_%s_%s' % (lane_num, fname))
            )

        fastq_paths = glob.glob(os.path.join(lane_dir, 'Undetermined_S0_L00%s_R*_001.fastq.gz' % lane_num))
        # move them up
        for fastq_path in fastq_paths:
            os.rename(fastq_path, os.path.join(fastq_dir, os.path.basename(fastq_path)))

    # then move the fastq associated with run elements
    for run_element in run_elements:
        fastq_paths = glob.glob(os.path.join(
            fastq_dir,
            'lane_%s' % run_element[ELEMENT_LANE],
            run_element[ELEMENT_PROJECT_ID],
            run_element[ELEMENT_SAMPLE_INTERNAL_ID],
            '*_L00%s_R*_001.fastq.gz' % run_element[ELEMENT_LANE]
        ))
        destination_dir = os.path.join(
            fastq_dir,
            run_element[ELEMENT_PROJECT_ID],
            run_element[ELEMENT_SAMPLE_INTERNAL_ID]
        )
        os.makedirs(destination_dir, exist_ok=True)
        for fastq_path in fastq_paths:
            os.rename(fastq_path, os.path.join(destination_dir, os.path.basename(fastq_path)))

    # All the files have been moved remove the lane directories
    for lane_dir in lane_dirs:
        shutil.rmtree(lane_dir)
