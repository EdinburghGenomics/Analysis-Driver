import os
import glob
import argparse
import sys

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from analysis_driver.quality_control.contamination_checks import ContaminationCheck
from analysis_driver.config import default as cfg


def main():
    args = _parse_args()
    args.func(args)

def _parse_args():
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers()

    sp_contamination_parser = subparsers.add_parser('contamination_check')
    sp_contamination_parser.add_argument('--fastq_files', required=True, nargs='+', help='the fastq file pairs')
    sp_contamination_parser.add_argument('--run_id', required=True)
    sp_contamination_parser.set_defaults(func=run_species_contamiantion_check)

    return parser.parse_args()

def run_species_contamiantion_check(args):
    species_contamination_check = ContaminationCheck(sorted(args.fastq_files),args.run_id)
    species_contamination_check.start()
    species_contamination_check.join()



if __name__ == '__main__':
    sys.exit(main())