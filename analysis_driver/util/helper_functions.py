# retrieving run_element_id, which differs between barcodeless and multiplexed runs
def get_run_element_id(run_id, lane_number, barcode=None):
    if barcode is None:
        return '{dataset_name}_{lane}'.format(dataset_name=run_id, lane=lane_number)
    else:
        return '{dataset_name}_{lane}_{barcode}'.format(dataset_name=run_id, lane=lane_number, barcode=barcode)
