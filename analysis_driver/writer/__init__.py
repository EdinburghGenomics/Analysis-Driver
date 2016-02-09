__author__ = 'mwham'
import csv
import os.path
from .script_writer import ScriptWriter
from .pbs_writer import PBSWriter
from . import bash_commands
from analysis_driver.app_logging import get_logger
from analysis_driver.config import default as cfg

app_logger = get_logger(__name__)


def write_bcbio_csv(run_dir, sample_id, fastqs, user_sample_id=None):
    """
    Write out a simple csv mapping fastq files to a sample id
    """
    if not user_sample_id:
        user_sample_id = sample_id
    csv_file = os.path.join(run_dir, 'samples_' + sample_id + '.csv')
    app_logger.info('Writing BCBio sample csv ' + csv_file)

    with open(csv_file, 'w') as f:
        writer = csv.writer(f)
        writer.writerow(['samplename', 'description'])
        for fq in fastqs:
            writer.writerow([fq, user_sample_id])

    return csv_file


def get_script_writer(job_name, run_id, walltime=None, cpus=1, mem=2, jobs=1, **kwargs):
    if cfg['job_execution'] == 'pbs':
        return PBSWriter(job_name, run_id, cpus, mem, walltime, jobs, **kwargs)
    # elif cfg['job_execution'] == 'sge':
        # return SGEWriter(job_name, run_id, walltime, cpus, mem, jobs)
    else:
        return ScriptWriter(job_name, run_id, jobs, **kwargs)
