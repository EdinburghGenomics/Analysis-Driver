import luigi
from os.path import join
from egcg_core import clarity
from egcg_core.app_logging import logging_default as log_cfg
from analysis_driver.dataset import RunDataset, SampleDataset
from analysis_driver.config import default as cfg
from analysis_driver.exceptions import PipelineError


class EGCGParameter(luigi.Parameter):
    """Parameter that does not call `warnings.warn` in self.serialize."""
    def serialize(self, x):
        return str(x)


class Stage(luigi.Task):
    __stagename__ = None
    exit_status = None
    previous_stages = ()
    dataset = EGCGParameter()

    @property
    def dataset_name(self):
        return self.dataset.name

    @property
    def job_dir(self):
        return join(cfg['jobs_dir'], self.dataset_name)

    @property
    def input_dir(self):
        return join(cfg.get('intermediate_dir', cfg['input_dir']), self.dataset_name)

    def get_cached_data(self, key):
        return self.dataset.data[key]

    def cache_data(self, key, value):
        self.dataset.data[key] = value

    def output(self):  # if <stage_name>.done is present, the stage is complete
        return RestAPITarget(self)

    def run(self):
        self.dataset.start_stage(self.stage_name)
        self.exit_status = self._run()
        self.dataset.end_stage(self.stage_name, self.exit_status)
        if self.exit_status:
            raise PipelineError('Exit status was %s. Stopping' % self.exit_status)

    @property
    def stage_name(self):
        if self.__stagename__:
            return self.__stagename__
        return self.__class__.__name__.lower()

    def _run(self):
        raise NotImplementedError

    def requires(self):
        """
        Generates prior Stages from self.previous_stages, which should be either a single luigi.Task or a
        tuple. If it's a tuple, each element can be a luigi.Task or a
        tuple[luigi.Task,dict[str,luigi.Parameter]].
        """
        if isinstance(self.previous_stages, type):
            return self.previous_stages(dataset=self.dataset)

        for s in self.previous_stages:
            if isinstance(s, type):
                yield s(dataset=self.dataset)
            elif type(s) is tuple:
                assert len(s) == 2 and type(s[1]) is dict, 'Invalid format for self.previous_stages'
                yield s[0](dataset=self.dataset, **s[1])


class RestAPITarget(luigi.Target):
    def __init__(self, stage):
        self.stage = stage

    @property
    def rest_api_stage(self):
        return self.stage.dataset.get_stage(self.stage.stage_name)

    def exists(self):
        s = self.rest_api_stage
        return s and bool(s.get('date_finished')) and s.get('exit_status') == 0


def pipeline(dataset):
    _setup_luigi_logging()

    if isinstance(dataset, RunDataset):
        from . import fastq_production
        final_stage = fastq_production.OutputQCData(dataset=dataset)
    # elif isinstance(dataset, SampleDataset):
    #     species = clarity.get_species_from_sample(dataset.name)
    #     if species is None:
    #         raise PipelineError('No species information found in the LIMS for ' + dataset.name)
    #     elif species == 'Homo sapiens':
    #         final_stage = BCBioVarCalling
    #     elif clarity.get_sample(dataset.name).udf.get('Analysis Type') == 'Variant Calling':
    #         final_stage = VarCalling
    #     else:
    #         final_stage = BasicQC
    # else:
    #     raise AssertionError('Unexpected dataset type: ' + str(dataset))

    luigi.build(
        [final_stage],
        local_scheduler=True
    )
    return final_stage.exit_status


def _setup_luigi_logging():
    luigi.interface.setup_interface_logging.has_run = True  # turn off Luigi's default logging setup
    log_cfg.get_logger('luigi-interface', 10)  # just calling log_cfg.get_logger registers the luigi-interface
