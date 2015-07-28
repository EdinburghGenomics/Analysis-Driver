__author__ = 'mwham'
import os
from analysis_driver.util import NamedAppLogger

app_logger = NamedAppLogger(__name__)


def find_fastqs(location, sample_project):
    app_logger.info('Looking for fastqs in ' + os.path.join(location, sample_project))

    fastq_dir = os.path.join(location, sample_project)
    fastqs = {}

    for sample_id in _list_dirs(fastq_dir):
        sample_id_path = os.path.join(fastq_dir, sample_id)

        fastqs[sample_id] = [
            os.path.join(sample_id_path, fq)
            for fq in os.listdir(sample_id_path)
            if fq.endswith('.fastq.gz')
        ]

    app_logger.info('Found ' + str(len(fastqs)) + ' fastq files')
    return fastqs


def flatten_fastqs(location, sample_projects):
    fastqs = []

    for sample_project in sample_projects:
        for fqs in find_fastqs(location, sample_project).values():
            fastqs = fastqs + fqs

    return fastqs


def _list_dirs(d):
    return [x for x in os.listdir(d) if os.path.isdir(os.path.join(d, x))]
