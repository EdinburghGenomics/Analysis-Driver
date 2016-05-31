import os
import csv
from glob import glob

import shutil

from analysis_driver import executor
from analysis_driver.exceptions import AnalysisDriverError, PipelineError
from analysis_driver.app_logging import logging_default as log_cfg
from analysis_driver.config import default as cfg
from analysis_driver.util.bash_commands import is_remote_path, rsync_from_to
from . import bash_commands

app_logger = log_cfg.get_logger('util')


def bcbio_prepare_samples_cmd(job_dir, sample_id, fastqs, user_sample_id=None):
    """
    Call bcbio_prepare_samples with a csv sample file and a list of fastqs.
    :param str job_dir: Full path to the run folder
    :param str sample_id: An ID to assign to the samples
    :param list fastqs: Full paths to each input fastq file
    """
    # setup the BCBio merged csv file
    bcbio_csv_file = _write_bcbio_csv(job_dir, sample_id, fastqs, user_sample_id=user_sample_id)
    app_logger.info('Setting up BCBio samples from ' + bcbio_csv_file)

    merged_dir = os.path.join(job_dir, 'merged')
    return ' '.join(
        (
            os.path.join(cfg['tools']['bcbio'], 'bin', 'bcbio_prepare_samples.py'),
            '--out',
            merged_dir,
            '--csv',
            bcbio_csv_file
        )
    )


def find_files(*path_parts):
    return glob(os.path.join(*path_parts))


def find_file(*path_parts):
    files = find_files(*path_parts)
    if files:
        return files[0]


def str_join(*parts, separator=''):
    return separator.join(parts)


def find_fastqs(location, sample_project, sample_id, lane=None):
    """
    Find all .fastq.gz files in an input folder 'location/sample_project'
    :param location: The overall directory to search
    :param sample_project: The sample_project directory to search
    :return: Full paths to *.fastq.gz files in the sample_project dir.
    :rtype: list[str]
    """
    if lane:
        pattern = os.path.join(sample_project, sample_id, '*L00%s*.fastq.gz' % lane)
    else:
        pattern = os.path.join(sample_project, sample_id, '*.fastq.gz')
    fastqs = find_files(location, pattern)
    app_logger.info('Found %s fastq files for %s', len(fastqs), pattern)
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
    app_logger.info('Found %s fastqs in %s', len(fastqs), location)
    return fastqs


def find_all_fastq_pairs(location):
    """
    Return the results of find_all_fastqs as a list or paired fastq files.
    It does not check the fastq name but expect that they will come together after sorting
    :return: Full paths to all *.fastq.gz files for all sample projects and sample ids in the input dir aggregated per pair
    :rtype: list[tuple(str, str)]
    """
    fastqs = find_all_fastqs(location)
    if len(fastqs) % 2 != 0 :
        raise AnalysisDriverError('Expect a even number of fastq file in %s found %s' % (location, len(fastqs)))
    fastqs.sort()
    return list(zip(*[iter(fastqs)]*2))




def _write_bcbio_csv(run_dir, sample_id, fastqs, user_sample_id=None):
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

def same_fs(file1, file2):
    if not file1 or not file2:
        return False
    if not os.path.exists(file1):
        return same_fs(os.path.dirname(file1), file2)
    if not os.path.exists(file2):
        return same_fs(file1, os.path.dirname(file2))

    dev1 = os.stat(file1).st_dev
    dev2 = os.stat(file2).st_dev
    return dev1 == dev2


def move_dir(src_dir, dest_dir):
    os.makedirs(dest_dir, exist_ok=True)
    list_file = os.listdir(src_dir)
    for f in list_file:
        if os.path.isdir(f):
            move_dir(os.path.join(src_dir, f), os.path.join(dest_dir, f))
        else:
            move_file_or_link_target(os.path.join(src_dir, f), dest_dir)
    return 0





def move_file_or_link_target(file_path, dest_dir):
    """move a file to a destination directory"""
    dest_file = os.path.join(dest_dir, os.path.basename(file_path))
    if os.path.islink(file_path):
        file_path = os.path.realpath(file_path)
    dest_file = shutil.move(file_path, dest_file)





