import os


def get_fastq_files(path):
    """
    Get the names of all fastq files in the given directory
    :param path
    :return: A list of fastq file names
    """
    # TODO: Check the paths exist
    assert os.path.isdir(path), 'Expected a valid directory to find fastq files'

    return [f for f in os.listdir(path) if f.endswith('.fastq')]
