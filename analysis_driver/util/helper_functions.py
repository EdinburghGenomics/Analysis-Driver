import os


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
