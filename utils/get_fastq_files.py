import os


def get_fastq_files(input_dir):
    """
    :param input_dir: Full path to an input directory
    :return: Names of all fastq files in input_dir
    """
    # TODO: Check the paths exist
    assert os.path.isdir(input_dir), 'Expected a valid directory to find fastq files'

    return [f for f in os.listdir(input_dir) if f.endswith('.fastq')]

