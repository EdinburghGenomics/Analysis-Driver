import os
import os.path

def make_project(path):
    """
    Create a temporary directory for a project
    :param path:
    :return: None
    """
    os.makedirs(path)
    os.makedirs(path + '/pbs')


def get_dir_name(path):
    return os.path.basename(os.path.normpath(path))
