import luigi
from os.path import join, exists
from egcg_core import clarity
from egcg_core.app_logging import AppLogger, logging_default as log_cfg
from analysis_driver.dataset import RunDataset, SampleDataset
from analysis_driver.config import default as cfg
from analysis_driver.exceptions import PipelineError

_dataset = None


class EGCGParameter(luigi.Parameter):
    """Parameter that does not call `warnings.warn` in self.serialize."""
    def serialize(self, x):
        if type(x) in (list, tuple, dict, set):
            x = '<%s of len %s>' % (x.__class__.__name__, len(x))
        return str(x)

    # def __getattribute__(self, item):
    #     return getattr(self, item, None)


class BasicStage(luigi.Task, AppLogger):
    __stagename__ = None
    previous_stages = ()

    @property
    def dataset(self):
        return _dataset

    @property
    def dataset_name(self):
        return self.dataset.name

    @property
    def job_dir(self):
        return self.dataset.job_dir

    @property
    def input_dir(self):
        return join(cfg.get('intermediate_dir', cfg['input_dir']), self.dataset_name)

    @property
    def stage_name(self):
        if self.__stagename__:
            return self.__stagename__
        return self.__class__.__name__.lower()

    def requires(self):
        """
        Generates prior Stages from self.previous_stages, which should be either a single luigi.Task or a
        tuple. If it's a tuple, each element can be a luigi.Task or a
        tuple[luigi.Task,dict[str,luigi.Parameter]].
        """
        p = self.previous_stages
        if type(p) is not tuple:
            p = (p,)

        for s in p:
            if isinstance(s, type):
                s = s()
            yield s


class Stage(BasicStage):
    exit_status = None
    expected_output_files = []

    def output(self):
        return [RestAPITarget(self)] + self.expected_output_files

    def run(self):
        self.dataset.start_stage(self.stage_name)
        self.exit_status = self._run()
        self.dataset.end_stage(self.stage_name, self.exit_status)
        if self.exit_status:
            raise PipelineError('Exit status was %s. Stopping' % self.exit_status)
        for f in self.expected_output_files:
            self.dataset.expected_output_files.append(f)
            self.debug('Registered output file %s' % f.filename)
        self.info('Finished with %s expected output files' % len(self.expected_output_files))

    def _run(self):
        raise NotImplementedError


class RestAPITarget(luigi.Target):
    def __init__(self, stage):
        self.stage = stage

    @property
    def rest_api_stage(self):
        return self.stage.dataset.get_stage(self.stage.stage_name)

    def exists(self):
        s = self.rest_api_stage
        return s and bool(s.get('date_finished')) and s.get('exit_status') == 0


class FileTarget(luigi.Target):
    def __init__(self, filename, deliver, new_basename=None, required=True):
        self.filename = filename
        self.new_basename = new_basename
        self.required = required
        self.deliver = deliver

    def exists(self):
        return exists(self.filename) or not self.required


def pipeline(dataset):
    _setup_luigi_logging()

    if isinstance(dataset, RunDataset):
        from .fastq_production import OutputQCData as FinalStage
    elif isinstance(dataset, SampleDataset):
        species = clarity.get_species_from_sample(dataset.name)
        if species is None:
            raise PipelineError('No species information found in the LIMS for ' + dataset.name)
        elif species == 'Homo sapiens':
            from .bcbio_var_calling import DataOutput as FinalStage
        elif clarity.get_sample(dataset.name).udf.get('Analysis Type') == 'Variant Calling':
            from .var_calling import VarCalling as FinalStage
        else:
            from .var_calling import BasicQC as FinalStage
    else:
        raise AssertionError('Unexpected dataset type: ' + str(dataset))

    global _dataset
    _dataset = dataset
    final_stage = FinalStage()
    luigi.build(
        [final_stage],
        local_scheduler=True
    )
    return final_stage.exit_status


def _setup_luigi_logging():
    luigi.interface.setup_interface_logging.has_run = True  # turn off Luigi's default logging setup
    log_cfg.get_logger('luigi-interface', 10)  # just calling log_cfg.get_logger registers the luigi-interface
