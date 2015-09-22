__author__ = 'mwham'
import os
from glob import glob
from analysis_driver.app_logging import get_logger

app_logger = get_logger('fastq_handler')


def find_fastqs(location, sample_project):
    """
    Iterate through an input folder and find all .fastq.gz files. The input folder should be
    'location/sample_project'
    :param location: The overall directory to search
    :param sample_project: The sample_project directory to search
    :return: A dict mapping sample ids to full paths to *.fastq.gz files in the sample_project dir.
    :rtype: dict[str, list[str]]
    """
    fastq_dir = os.path.join(location, sample_project)
    fastqs = {}

    for sample_id in os.listdir(fastq_dir):
        fastqs[sample_id] = glob(os.path.join(fastq_dir, sample_id, '*.fastq.gz'))

    app_logger.info('Found %s fastq files in %s' % (len(fastqs), os.path.join(location, sample_project)))
    return fastqs


def flatten_fastqs(location):
    """
    Return the results of find_fastqs as a flat list.
    :return: Full paths to all *.fastq.gz files for all sample projects and sample ids in the input dir
    :rtype: list[str]
    """
    app_logger.debug('Flattening fastqs')
    return glob(os.path.join(location, '*', '*', '*.fastq.gz'))
