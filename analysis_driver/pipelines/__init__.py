import luigi
from egcg_core import clarity
from egcg_core.app_logging import logging_default as log_cfg
from egcg_core.config import cfg
from analysis_driver import segmentation
from analysis_driver.dataset_scanner import RunDataset, SampleDataset, ProjectDataset
from analysis_driver.exceptions import PipelineError
from analysis_driver.pipelines.qc_pipelines import qc_pipeline
from analysis_driver.pipelines.bcbio_pipelines import bcbio_var_calling_pipeline
from analysis_driver.pipelines.demultiplexing import demultiplexing_pipeline
from analysis_driver.pipelines.variant_calling import var_calling_pipeline
from analysis_driver.pipelines.projects import project_pipeline


def pipeline(d):
    _setup_luigi_logging()

    if isinstance(d, RunDataset):
        final_stage = Demultiplexing(dataset=d)
    elif isinstance(d, SampleDataset):
        species = clarity.get_species_from_sample(d.name)
        analysis_type = clarity.get_sample(d.name).udf.get('Analysis Type')
        if species is None:
            raise PipelineError('No species information found in the LIMS for ' + d.name)
        elif species == 'Homo sapiens':
            final_stage = BCBioVarCalling(dataset=d, species=species, analysis_type=analysis_type)
        elif analysis_type == 'Variant Calling':
            final_stage = VarCalling(dataset=d, species=species)
        else:
            final_stage = QC(dataset=d, species=species)
    elif isinstance(d, ProjectDataset):
        final_stage = Project(dataset=d)
    else:
        raise AssertionError('Unexpected dataset type: ' + str(d))

    final_stage.exit_status = 9

    luigi_params = {'tasks': [final_stage], 'local_scheduler': cfg.query('luigi', 'local_scheduler')}
    if luigi_params['local_scheduler'] is not True:
        luigi_params['scheduler_url'] = cfg['luigi']['scheduler_url']

    luigi.build(**luigi_params)
    return final_stage.exit_status


def _setup_luigi_logging():
    luigi.interface.setup_interface_logging.has_run = True  # turn off Luigi's default logging setup
    log_cfg.get_logger('luigi-interface', 10)  # just calling log_cfg.get_logger registers the luigi-interface


class Demultiplexing(segmentation.BasicStage):
    def run(self):
        self.exit_status = demultiplexing_pipeline(self.dataset)


class BCBioVarCalling(segmentation.BasicStage):
    species = segmentation.EGCGParameter()
    analysis_type = segmentation.EGCGParameter()

    def run(self):
        self.exit_status = bcbio_var_calling_pipeline(self.dataset, self.species, self.analysis_type)


class VarCalling(segmentation.BasicStage):
    species = segmentation.EGCGParameter()

    def run(self):
        self.exit_status = var_calling_pipeline(self.dataset, self.species)


class QC(segmentation.BasicStage):
    species = segmentation.EGCGParameter()

    def run(self):
        self.exit_status = qc_pipeline(self.dataset, self.species)


class Project(segmentation.BasicStage):
    def run(self):
        self.exit_status = project_pipeline(self.dataset)
