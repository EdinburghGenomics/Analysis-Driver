__author__ = 'mwham'
import os.path
from analysis_driver.config import default as cfg
from .script_writer import ScriptWriter
from .pbs_writer import PBSWriter


def get_script_writer(job_name, run_id, walltime=None, cpus=1, mem=2, jobs=1, **kwargs):
    working_dir = os.path.join(cfg['jobs_dir'], run_id)
    if cfg['job_execution'] == 'pbs':
        return PBSWriter(job_name, working_dir, cfg['job_queue'], cpus, mem, walltime, jobs, **kwargs)
    # elif cfg['job_execution'] == 'sge':
        # return SGEWriter(job_name, working_dir, cfg['job_queue'], walltime, cpus, mem, jobs)
    else:
        return ScriptWriter(job_name, working_dir, cfg['job_queue'], jobs, **kwargs)
