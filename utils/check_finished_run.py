import os.path

def file_exists(file):
    """
    Check if a file has been created
    :param file: Full path to the file
    :return: True if exists, False if not created yet
    """
    return os.path.exists(file)


# Unit Test
if __name__ == '__main__':
    filename = '/home/U008/lcebaman/scripts/data/Finish.txt'
    print(file_exists(filename))
