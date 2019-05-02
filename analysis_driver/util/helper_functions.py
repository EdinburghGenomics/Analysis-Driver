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
    It recurse if one of the value is found to be a dict.
    This function does not modify the provided data_structure.
    """
    new_data_structure = {}
    for k in data_structure:
        if isinstance(data_structure[k], dict):
            new_data_structure[k] = prepend_path_to_data_files(prepend_path, data_structure[k])
        else:
            new_data_structure[k] = os.path.join(prepend_path, data_structure[k])
    return new_data_structure
