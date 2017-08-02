import luigi
from egcg_core import clarity
from egcg_core.app_logging import logging_default as log_cfg
from egcg_core.config import cfg
from analysis_driver.dataset_scanner import RunDataset, SampleDataset, ProjectDataset
from analysis_driver.exceptions import PipelineError
from analysis_driver.pipelines import demultiplexing, bcbio, qc, variant_calling, projects
from analysis_driver.tool_versioning import toolset


def pipeline(dataset):
    """
    Decide which pipeline to run on a dataset, and run luigi.build.
    :param analysis_driver.dataset.Dataset dataset:
    """
    luigi.interface.setup_interface_logging.has_run = True  # turn off Luigi's default logging setup
    log_cfg.get_logger('luigi-interface', 20)  # just calling log_cfg.get_logger registers the luigi-interface

    if isinstance(dataset, RunDataset):
        _pipeline = demultiplexing
        _pipeline_type = 'run_processing'

    elif isinstance(dataset, SampleDataset):
        analysis_type = clarity.get_sample(dataset.name).udf.get('Analysis Type')
        if dataset.species is None:
            raise PipelineError('No species information found in the LIMS for ' + dataset.name)

        elif dataset.species == 'Homo sapiens':
            _pipeline_type = 'human_sample_processing'
            _pipeline = bcbio
        elif analysis_type in ['Variant Calling', 'Variant Calling gatk']:
            _pipeline_type = 'non_human_sample_processing'
            _pipeline = variant_calling
        else:
            _pipeline_type = 'non_human_sample_processing'
            _pipeline = qc
    elif isinstance(dataset, ProjectDataset):
        _pipeline_type = 'project_processing'
        _pipeline = projects
    else:
        raise AssertionError('Unexpected dataset type: ' + str(dataset))

    dataset.pipeline_type = _pipeline_type
    toolset.select_type(_pipeline_type)
    toolset.select_version(dataset.toolset_version)

    dataset.start()
    final_stage = _pipeline.build_pipeline(dataset)

    luigi_params = {
        'tasks': [final_stage],
        'local_scheduler': cfg.query('luigi', 'local_scheduler'),
        'workers': cfg.query('luigi', 'max_parallel_jobs', ret_default=4)
    }
    if luigi_params['local_scheduler'] is not True:
        luigi_params['scheduler_url'] = cfg['luigi']['scheduler_url']

    success = luigi.build(**luigi_params)

    # if any exception occurred during the pipeline raise them here again
    dataset.raise_exceptions()

    return 0 if success is True else 9
