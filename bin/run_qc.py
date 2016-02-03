__author__ = 'tcezard'
import sys
import os
import argparse
import logging

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from analysis_driver.quality_control.genotype_validation import GenotypeValidation
from analysis_driver.config import logging_default as log_cfg
log_cfg.default_level = logging.DEBUG
log_cfg.add_handler('stdout', logging.StreamHandler(stream=sys.stdout), logging.DEBUG)


def main():
    args = _parse_args()
    args.func(args)

def _parse_args():
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers()

    run_parser = subparsers.add_parser('genotype')
    run_parser.add_argument('--fastqs_files', nargs='+')
    run_parser.add_argument('--sample_id', required = True)
    run_parser.set_defaults(func=run_genotype_validation)

    return parser.parse_args()


def run_genotype_validation(args):
    geno_val = GenotypeValidation(sorted(args.fastqs_files), args.sample_id)
    geno_val.start()
    validation_results = geno_val.join()

if __name__ == '__main__':
    sys.exit(main())
