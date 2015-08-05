__author__ = 'mwham'

import sys
import os.path
import argparse

if __name__ == '__main__':

    sys.path.append(os.path.dirname(os.path.dirname(__file__)))
    import analysis_driver

    parser = argparse.ArgumentParser()
    parser.add_argument('input_run_folder', type=str, help='An absolute path to an input data directory')
    # TODO: add a --log-level flag

    args = parser.parse_args(sys.argv[1:])

    sys.exit(analysis_driver.driver.main(input_run_folder=args.input_run_folder))
