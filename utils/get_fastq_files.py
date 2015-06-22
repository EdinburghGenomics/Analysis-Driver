import os

def get_fastq_files(path):
    """
    Return a list with the names of fastq files within the given directory
    :param path:
    :return:
    """
    # TODO: Check the path exists
    assert os.path.isdir(path), 'Expected a valid directory to find fastq files'

    return [file for file in os.listdir(path) if file.endswith('.fastq')]
