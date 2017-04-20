import luigi
from egcg_core import clarity
from egcg_core.app_logging import logging_default as log_cfg
from egcg_core.config import cfg
from analysis_driver.dataset_scanner import RunDataset, SampleDataset, ProjectDataset
from analysis_driver.exceptions import PipelineError
from analysis_driver.pipelines import demultiplexing, bcbio, qc, variant_calling, projects


def pipeline(d):
    luigi.interface.setup_interface_logging.has_run = True  # turn off Luigi's default logging setup
    log_cfg.get_logger('luigi-interface', 10)  # just calling log_cfg.get_logger registers the luigi-interface

    if isinstance(d, RunDataset):
        _pipeline = demultiplexing
    elif isinstance(d, SampleDataset):
        species = clarity.get_species_from_sample(d.name)
        analysis_type = clarity.get_sample(d.name).udf.get('Analysis Type')
        if species is None:
            raise PipelineError('No species information found in the LIMS for ' + d.name)
        elif species == 'Homo sapiens':
            _pipeline = bcbio
        elif analysis_type == 'Variant Calling':
            _pipeline = variant_calling
        else:
            _pipeline = qc
    elif isinstance(d, ProjectDataset):
        _pipeline = projects
    else:
        raise AssertionError('Unexpected dataset type: ' + str(d))

    final_stage = _pipeline.build_pipeline(d)
    final_stage.exit_status = 9

    luigi_params = {'tasks': [final_stage], 'local_scheduler': cfg.query('luigi', 'local_scheduler')}
    if luigi_params['local_scheduler'] is not True:
        luigi_params['scheduler_url'] = cfg['luigi']['scheduler_url']

    luigi.build(**luigi_params)
    return final_stage.exit_status
