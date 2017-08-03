import luigi
from egcg_core.config import cfg
from egcg_core.app_logging import logging_default as log_cfg
from analysis_driver.pipelines import demultiplexing, bcbio, qc, variant_calling, projects


def pipeline(dataset):
    """
    Decide which pipeline to run on a dataset, and run luigi.build.
    :param analysis_driver.dataset.Dataset dataset:
    """
    luigi.interface.setup_interface_logging.has_run = True  # turn off Luigi's default logging setup
    log_cfg.get_logger('luigi-interface', 20)  # just calling log_cfg.get_logger registers the luigi-interface

    dataset.resolve_pipeline_and_toolset()
    dataset.start()
    final_stage = dataset.pipeline.build_pipeline(dataset)

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
