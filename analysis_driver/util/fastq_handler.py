__author__ = 'mwham'
import os
from glob import glob
from analysis_driver.app_logging import get_logger

app_logger = get_logger('fastq_handler')


def find_fastqs(location, sample_project, sample_id, lane=None):
    """
    Iterate through an input folder and find all .fastq.gz files. The input folder should be
    'location/sample_project'
    :param location: The overall directory to search
    :param sample_project: The sample_project directory to search
    :return: A list of full paths to *.fastq.gz files in the sample_project dir.
    :rtype: list[str]
    """
    if lane:
        pattern = os.path.join(sample_project, sample_id, '*L00%s*.fastq.gz' % lane)
    else:
        pattern = os.path.join(sample_project, sample_id, '*.fastq.gz')
    fastqs = glob(os.path.join(location, pattern))
    app_logger.info('Found %s fastq files for %s' % (len(fastqs), pattern))
    return fastqs


def find_all_fastqs(location):
    """
    Return the results of find_fastqs as a flat list.
    :return: Full paths to all *.fastq.gz files for all sample projects and sample ids in the input dir
    :rtype: list[str]
    """
    fastqs = []
    for name, dirs, files in os.walk(location):
        fastqs.extend([os.path.join(name, f) for f in files if f.endswith('.fastq.gz')])
    app_logger.info('Found %s fastqs in %s' % (len(fastqs), location))
    return fastqs
