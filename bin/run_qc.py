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
from analysis_driver.quality_control.contamination_checks import ContaminationCheck
from analysis_driver.reader.demultiplexing_parsers import parse_fastqscreen_file
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

    sp_contamination_parser = subparsers.add_parser('contamination_check')
    sp_contamination_parser.add_argument('--fastq_files', required=True, nargs='+', help='the fastq file pairs')
    sp_contamination_parser.add_argument('--run_id', required=True)
    sp_contamination_parser.add_argument('--sample_id', required=True)
    sp_contamination_parser.set_defaults(func=run_species_contamiantion_check)

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
    seq_vcf_file, validation_results = geno_val.join()
    user_sample_id = get_user_sample_name(sample_name=args.sample_id)
    output_commands = []
    for f in [seq_vcf_file, validation_results]:
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


def run_species_contamiantion_check(args):
    work_dir = os.path.join(cfg['jobs_dir'], args.run_id)
    os.makedirs(work_dir,exist_ok=True)
    species_contamination_check = ContaminationCheck(sorted(args.fastq_files),args.run_id)
    species_contamination_check.start()
    expected_output_files = species_contamination_check.join()
    expected_output_files = (''.join(expected_output_files))
    fastqscreen_result = parse_fastqscreen_file(expected_output_files, args.sample_id)
    print(fastqscreen_result)





if __name__ == '__main__':
    sys.exit(main())
