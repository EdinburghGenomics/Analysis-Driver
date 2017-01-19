import luigi
from egcg_core import clarity
from egcg_core.app_logging import logging_default as log_cfg
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
        final_stage = Demultiplexing
    elif isinstance(d, SampleDataset):
        species = clarity.get_species_from_sample(d.name)
        if species is None:
            raise PipelineError('No species information found in the LIMS for ' + d.name)
        elif species == 'Homo sapiens':
            final_stage = BCBioVarCalling
        elif clarity.get_sample(d.name).udf.get('Analysis Type') == 'Variant Calling':
            final_stage = VarCalling
        else:
            final_stage = QC
    elif isinstance(d, ProjectDataset):
        final_stage = Project
    else:
        raise AssertionError('Unexpected dataset type: ' + str(d))

    segmentation.dataset = d
    final_stage = final_stage()
    luigi.build(
        [final_stage],
        local_scheduler=True
    )
    return final_stage.exit_status


def _setup_luigi_logging():
    luigi.interface.setup_interface_logging.has_run = True  # turn off Luigi's default logging setup
    log_cfg.get_logger('luigi-interface', 10)  # just calling log_cfg.get_logger registers the luigi-interface


class Demultiplexing(segmentation.BasicStage):
    def run(self):
        self.exit_status = demultiplexing_pipeline(self.dataset)


class BCBioVarCalling(segmentation.BasicStage):
    def run(self):
        s = clarity.get_sample(self.dataset.name)
        self.exit_status = bcbio_var_calling_pipeline(
            self.dataset, s.udf.get('Genome Version'), s.udf.get('Analysis Type')
        )


class VarCalling(segmentation.BasicStage):
    def run(self):
        self.exit_status = var_calling_pipeline(
            self.dataset, clarity.get_species_from_sample(self.dataset.name)
        )


class QC(segmentation.BasicStage):
    def run(self):
        self.exit_status = qc_pipeline(self.dataset, clarity.get_species_from_sample(self.dataset.name))


class Project(segmentation.BasicStage):
    def run(self):
        self.exit_status = project_pipeline(self.dataset)
