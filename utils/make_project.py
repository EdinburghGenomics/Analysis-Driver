import os
from os.path import normpath, basename


def make_project(input_path):
    """
    Create a temporary directory for each PBS project

    :param input_path:
    :return:
    """
    os.makedirs(input_path)
    os.makedirs(input_path + '/pbs')


def get_dirname(input_path):
    return basename(normpath(input_path))
