__author__ = 'mwham'
from analysis_driver.config import default as cfg
from .script_writer import ScriptWriter
from .pbs_writer import PBSWriter


def get_script_writer(job_name, working_dir, walltime=None, cpus=1, mem=2, jobs=1, **kwargs):
    if cfg['job_execution'] == 'pbs':
        return PBSWriter(job_name, working_dir, cfg['job_queue'], cpus, mem, walltime, jobs, **kwargs)
    # elif cfg['job_execution'] == 'sge':
        # return SGEWriter(job_name, working_dir, cfg['job_queue'], walltime, cpus, mem, jobs)
    else:
        return ScriptWriter(job_name, working_dir, cfg['job_queue'], jobs, **kwargs)
