from egcg_core import clarity
from analysis_driver.dataset_scanner import RunDataset, SampleDataset
from analysis_driver.exceptions import PipelineError
from egcg_core.app_logging import logging_default as log_cfg

from analysis_driver.pipelines.qc_pipelines import qc_pipeline
from analysis_driver.pipelines.bcbio_pipelines import bcbio_var_calling_pipeline
from analysis_driver.pipelines.demultiplexing import demultiplexing_pipeline
from analysis_driver.pipelines.variant_calling import var_calling_pipeline

app_logger = log_cfg.get_logger('pipelines')


def pipeline(dataset):

    if isinstance(dataset, RunDataset):
        return demultiplexing_pipeline(dataset)
    elif isinstance(dataset, SampleDataset):
        species = clarity.get_species_from_sample(dataset.name)
        genome_version = clarity.get_sample(dataset.name).udf.get('Genome Version')
        analysis_type = clarity.get_sample(dataset.name).udf.get('Analysis Type')
        if species is None:
            raise PipelineError('No species information found in the LIMS for ' + dataset.name)
        elif species == 'Homo sapiens':
            return bcbio_var_calling_pipeline(dataset, genome_version, analysis_type)
        elif clarity.get_sample(dataset.name).udf.get('Analysis Type') == 'Variant Calling':
            return var_calling_pipeline(dataset, species)
        else:
            return qc_pipeline(dataset, species)
    else:
        raise AssertionError('Unexpected dataset type: ' + str(dataset))

