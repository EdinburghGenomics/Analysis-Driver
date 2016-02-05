import glob
import sys
import os
import argparse
import logging

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from analysis_driver.quality_control.genotype_validation import GenotypeValidation
from analysis_driver.config import logging_default as log_cfg
from analysis_driver import executor
from analysis_driver.exceptions import AnalysisDriverError
from analysis_driver.writer.bash_commands import rsync_from_to, is_remote_path
from analysis_driver.config import default as cfg

log_cfg.default_level = logging.DEBUG
log_cfg.add_handler('stdout', logging.StreamHandler(stream=sys.stdout), logging.DEBUG)


def main():
    args = _parse_args()
    args.func(args)

def _parse_args():
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers()

    geno_val_parser = subparsers.add_parser('genotype_validation')
    geno_val_parser.add_argument('--project_id', required=True)
    geno_val_parser.add_argument('--sample_id', required = True)
    geno_val_parser.set_defaults(func=run_genotype_validation)

    return parser.parse_args()


def run_genotype_validation(args):
    #Get the sample specific config
    cfg.merge(cfg['sample'])
    projects_source = cfg.query('output_dir')
    work_dir = os.path.join(cfg['jobs_dir'], args.sample_id)
    os.makedirs(work_dir,exist_ok=True)
    #Hack to retrive the fastq file from the CIFS share
    if is_remote_path(projects_source):
        # Need to retrieve the data localy
        fastq_files = os.path.join(projects_source, args.project_id, args.sample_id, '*_R?.fastq.gz')
        cmd = rsync_from_to(fastq_files, work_dir)
        exit_status = executor.execute(
                [cmd],
                job_name='getfastq',
                run_id=args.sample_id,
                cpus=1,
                mem=2
        ).join()
        if exit_status != 0:
            raise AnalysisDriverError("Copy of the fastq files from remote has failed")

        fastq_files = glob.glob(os.path.join(work_dir, '*_R?.fastq.gz'))
    else:
        fastq_files = glob.glob(os.path.join(projects_source, args.project_id, args.sample_id, '*_R?.fastq.gz'))

    geno_val = GenotypeValidation(sorted(fastq_files), args.sample_id)
    geno_val.start()
    validation_results = geno_val.join()

if __name__ == '__main__':
    sys.exit(main())
