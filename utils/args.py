import sys
import getopt
import os

def parse_args(argv):

    input_file = ''
 
    try:
        opts, args = getopt.getopt(argv, 'hi:', ['ifile='])
    except getopt.GetoptError:
        print('test.py -i <input_file> ')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print('test.py -i <input_file> ')
            sys.exit()
        elif opt in ('-i', '--ifile'):
            input_file = arg
            assert os.path.exists(input_file), 'Non existing directory, ' + input_file
        else:
            print('test.py -i <input_file> ')
            sys.exit()

    # print('Input directory is', input_file.rstrip('/')
    return input_file.rstrip('/')