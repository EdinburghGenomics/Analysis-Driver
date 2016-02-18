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
from analysis_driver.clarity import get_user_sample_name

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
    geno_val_parser.add_argument('--check_plate', action='store_true', default=False)
    geno_val_parser.add_argument('--check_project', action='store_true', default=False)
    geno_val_parser.set_defaults(func=run_genotype_validation)

    return parser.parse_args()


def run_genotype_validation(args):

    def retrive_data(paths, work_dir, sample_id, allow_fail=False):
        cmd = rsync_from_to(paths, work_dir)
        exit_status = executor.execute(
                [cmd],
                job_name='retrieve_data',
                run_id=sample_id,
                cpus=1,
                mem=2
        ).join()
        if exit_status != 0 and not allow_fail:
            raise AnalysisDriverError("Copy of the file(s) from remote has failed")

    #Get the sample specific config
    cfg.merge(cfg['sample'])
    projects_source = cfg.query('output_dir')
    work_dir = os.path.join(cfg['jobs_dir'], args.sample_id)
    os.makedirs(work_dir,exist_ok=True)
    #Hack to retrive the fastq file from the CIFS share
    if is_remote_path(projects_source):
        # First try to retrieve the genotype vcf file
        genotype_vcfs = os.path.join(projects_source, args.project_id, args.sample_id, '*_genotype_validation.vcf.gz')
        retrive_data(genotype_vcfs, work_dir, args.sample_id, allow_fail=True)
        genotype_vcfs = glob.glob(os.path.join(work_dir, '*_genotype_validation.vcf.gz'))
        if not genotype_vcfs:
            # Need to retrieve the fastq files localy
            fastq_files = os.path.join(projects_source, args.project_id, args.sample_id, '*_R?.fastq.gz')
            retrive_data(fastq_files, work_dir, args.sample_id)
            fastq_files = glob.glob(os.path.join(work_dir, '*_R?.fastq.gz'))
            genotype_vcf = None
        else:
            genotype_vcf = genotype_vcfs[0]
            fastq_files = []
    else:
        genotype_vcfs = glob.glob(os.path.join(work_dir, '*_genotype_validation.vcf.gz'))
        if not genotype_vcfs:
            fastq_files = glob.glob(os.path.join(projects_source, args.project_id, args.sample_id, '*_R?.fastq.gz'))
            genotype_vcf = None
        else:
            genotype_vcf = genotype_vcfs[0]
            fastq_files = []

    geno_val = GenotypeValidation(sorted(fastq_files), args.sample_id, vcf_file=genotype_vcf,
                                  check_plate=args.check_plate, check_project=args.check_project)
    geno_val.start()

    seq_vcf_file, sample2validation_results = geno_val.join()
    user_sample_id = get_user_sample_name(sample_name=args.sample_id)
    output_commands = []
    for f in [seq_vcf_file, sample2validation_results.get(args.sample_id)]:
        if f and os.path.exists(f):
            bf = os.path.basename(f)
            out_file = os.path.join(projects_source, args.project_id, args.sample_id, bf.replace(args.sample_id, user_sample_id))
            output_commands.append(rsync_from_to(f, out_file))

    exit_status = executor.execute(
            output_commands,
            job_name='output_results',
            run_id=args.sample_id,
            cpus=1,mem=2).join()

    if exit_status != 0:
        raise AnalysisDriverError("Copy of the results files to remote has failed")


if __name__ == '__main__':
    sys.exit(main())
