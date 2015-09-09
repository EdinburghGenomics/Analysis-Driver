__author__ = 'mwham'
from .script_writer import ScriptWriter
from .pbs_writer import PBSWriter
from . import commands
import csv
import os.path
from analysis_driver.app_logging import get_logger
from analysis_driver.config import default as cfg

app_logger = get_logger(__name__)


def write_bcbio_csv(run_dir, sample_id, fastqs):
    csv_file = os.path.join(run_dir, 'samples_' + sample_id + '.csv')
    app_logger.info('Writing BCBio sample csv ' + csv_file)

    with open(csv_file, 'w') as f:
        writer = csv.writer(f)

        writer.writerow(['samplename', 'description'])
        for fq in fastqs:
            writer.writerow([fq, sample_id])

    return csv_file


def get_script_writer(job_name, run_id, walltime, cpus, mem, jobs=1):
    if cfg['job_execution'] == 'pbs':
        return PBSWriter(job_name, run_id, walltime, cpus, mem, jobs)
    else:
        return ScriptWriter(job_name, run_id, jobs)


def write_jobs(writer, jobs, log_file_base=None):
    """
    :param ScriptWriter writer:
    :param list jobs: a list of commands to be written
    :return: the name of the script written
    """
    if len(jobs) == 1:
        writer.write_line(jobs[0])
    else:
        writer.start_array()
        for idx, job in enumerate(jobs):
            if log_file_base:
                writer.write_array_cmd(idx + 1, job, log_file_base + str(idx + 1) + '.log')
            else:
                writer.write_array_cmd(idx + 1, job)
        writer.finish_array()
    writer.save()

    return writer.script_name
