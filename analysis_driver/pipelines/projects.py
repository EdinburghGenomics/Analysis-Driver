from egcg_core import rest_communication


def project_pipeline(dataset):
    project_id = dataset.name
    project = rest_communication.get_document('projects', where={'project_id': project_id})
    # run the project pipeline
    exit_status = 0
    return exit_status
