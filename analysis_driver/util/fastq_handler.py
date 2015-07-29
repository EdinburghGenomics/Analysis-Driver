__author__ = 'mwham'
import os
from analysis_driver.util import NamedAppLogger

app_logger = NamedAppLogger(__name__)


def find_fastqs(location, sample_project):
    """
    Iterate through an input folder and find all .fastq.gz files.
    The function drills down directly to a sample project. Run the function in a second layer of iteration to
    get all fastqs in a fastqc output folder.
    :param location: The overall directory to search
    :param sample_project: The sample_project directory to search
    :return: A dict mapping sample ids to full paths to *.fastq.gz files in the sample_project dir.
    :rtype: dict[str, list[str]]
    """
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
    """
    Find all *.fastq.gz files as a flat list. Calls find_fastqs as part of its functionality.
    :param str location: A file path to the input_data folder containing
    :param list/tuple sample_projects: A list of sample projects to search in the fastq dir.
    :return: Full paths to all *.fastq.gz files in the input dir
    :rtype: list[str]
    """
    fastqs = []

    for sample_project in sample_projects:
        for fqs in find_fastqs(location, sample_project).values():
            fastqs = fastqs + fqs

    return fastqs


def _list_dirs(d):
    """
    Find all subdirectories in an input dir.
    :param d: A directory to search
    :return: A list of all subdirectories (not full paths)
    :rtype: list[str]
    """
    return [x for x in os.listdir(d) if os.path.isdir(os.path.join(d, x))]
