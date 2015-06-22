import os.path
import sys

def file_exists(file):
    """
    Check that a file has been created
    :arg file: A path to the input file
    :return: True (file exists) or False (file has not been created yet)
    """
    try:
        exists = os.path.isfile(file) 
    except os.error:
        print('ERROR in fileExists')
        sys.exit()
   
    return exists


if __name__ == '__main__':
    # Unit Test
    filename = "/home/U008/lcebaman/scripts/data/Finish.txt"
    print(file_exists(filename))
