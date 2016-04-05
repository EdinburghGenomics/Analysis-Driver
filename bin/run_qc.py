import glob
import sys
import os
import argparse
import logging

from analysis_driver.dataset_scanner import SampleScanner
from analysis_driver.quality_control.gender_validation import GenderValidation

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from analysis_driver.quality_control.genotype_validation import GenotypeValidation
from analysis_driver.app_logging import logging_default as log_cfg
from analysis_driver import executor
from analysis_driver.exceptions import AnalysisDriverError
from analysis_driver.util.bash_commands import rsync_from_to, is_remote_path
from analysis_driver.config import default as cfg
from analysis_driver.clarity import get_user_sample_name
from analysis_driver.quality_control.contamination_checks import ContaminationCheck, VerifyBamId
from analysis_driver.reader.demultiplexing_parsers import get_fastqscreen_results
log_cfg.default_level = logging.DEBUG
log_cfg.add_handler(logging.StreamHandler(stream=sys.stdout), logging.DEBUG)


def main():
    args = _parse_args()
    args.func(args)


def _parse_args():
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers()

    geno_val_parser = subparsers.add_parser('genotype_validation')
    geno_val_parser.add_argument('--project_id', required=True)
    geno_val_parser.add_argument('--sample_id', required = True)
    geno_val_parser.add_argument('--check_neighbour', action='store_true', default=False)
    geno_val_parser.add_argument('--check_project', action='store_true', default=False)
    geno_val_parser.add_argument('--check_samples', nargs='*')
    geno_val_parser.set_defaults(func=run_genotype_validation)

    sp_contamination_parser = subparsers.add_parser('species_contamination_check')
    sp_contamination_parser.add_argument('--fastq_files', required=True, nargs='+', help='the fastq file pairs')
    sp_contamination_parser.add_argument('--work_dir', required=True)
    sp_contamination_parser.add_argument('--sample_id', required=True)
    sp_contamination_parser.set_defaults(func=run_species_contamination_check)

    sample_contamination_parser = subparsers.add_parser('sample_contamination_check')
    sample_contamination_parser.add_argument('--bam_file', required=True, help='the fastq file pairs')
    sample_contamination_parser.add_argument('--work_dir', required=False)
    sample_contamination_parser.add_argument('--sample_id', required=True)
    sample_contamination_parser.set_defaults(func=run_sample_contamination_check)

    gender_valid_parser = subparsers.add_parser('gender_validation')
    gender_valid_parser.add_argument('--sample_id', type=str, help='sample ID for creating a Sample dataset object')
    gender_valid_parser.add_argument('-v', '--vcf_file', dest="vcf_file", type=str, help='the vcf file used to detect the gender')
    gender_valid_parser.add_argument('-s', '--working_dir', dest="working_dir", type=str, help='the working dir for execution')

    return parser.parse_args()


def run_genotype_validation(args):
    def retrieve_data(paths, work_dir, allow_fail=False):
        cmd = rsync_from_to(paths, work_dir)
        exit_status = executor.execute(
                [cmd],
                job_name='retrieve_data',
                working_dir=work_dir,
                cpus=1,
                mem=2
        ).join()
        if exit_status != 0 and not allow_fail:
            raise AnalysisDriverError("Copy of the file(s) from remote has failed")

    # Get the sample specific config
    cfg.merge(cfg['sample'])
    projects_source = cfg.query('output_dir')
    work_dir = os.path.join(cfg['jobs_dir'], args.sample_id)
    os.makedirs(work_dir, exist_ok=True)
    dataset = SampleScanner(cfg).get_dataset(args.sample_id)


    # Hack to retrieve the fastq file from the CIFS share
    if is_remote_path(projects_source):
        # First try to retrieve the genotype vcf file
        genotype_vcfs = os.path.join(projects_source, args.project_id, args.sample_id, '*_genotype_validation.vcf.gz')
        retrieve_data(genotype_vcfs, work_dir, allow_fail=True)
        genotype_vcfs = glob.glob(os.path.join(work_dir, '*_genotype_validation.vcf.gz'))
        if not genotype_vcfs:
            # Need to retrieve the fastq files localy
            fastq_files = os.path.join(projects_source, args.project_id, args.sample_id, '*_R?.fastq.gz')
            retrieve_data(fastq_files, work_dir)
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

    geno_val = GenotypeValidation(dataset, work_dir, sorted(fastq_files), args.sample_id, vcf_file=genotype_vcf,
                                  check_neighbour=args.check_neighbour, check_project=args.check_project,
                                  list_samples=args.check_samples)
    geno_val.start()

    seq_vcf_file, validation_results = geno_val.join()
    user_sample_id = get_user_sample_name(sample_name=args.sample_id)
    output_commands = []
    for f in [seq_vcf_file, validation_results]:
        if f and os.path.exists(f):
            bf = os.path.basename(f)
            out_file = os.path.join(projects_source, args.project_id, args.sample_id, bf.replace(args.sample_id, user_sample_id))
            output_commands.append(rsync_from_to(f, out_file))

    exit_status = executor.execute(
        output_commands,
        job_name='output_results',
        working_dir=work_dir,
        cpus=1,
        mem=2
    ).join()

    if exit_status != 0:
        raise AnalysisDriverError("Copy of the results files to remote has failed")


def run_species_contamination_check(args):
    if args.work_dir:
        work_dir = args.work_dir
    else:
        work_dir = os.path.join(cfg['jobs_dir'], args.sample_id)
    os.makedirs(work_dir, exist_ok=True)
    dataset = SampleScanner(cfg).get_dataset(args.sample_id)
    species_contamination_check = ContaminationCheck(dataset, work_dir, sorted(args.fastq_files))
    species_contamination_check.start()
    expected_output_files = species_contamination_check.join()
    expected_output_files = (''.join(expected_output_files))
    fastqscreen_result = get_fastqscreen_results(expected_output_files, args.sample_id)
    print(fastqscreen_result)

def run_sample_contamination_check(args):
    if args.work_dir:
        work_dir = args.work_dir
    else:
        work_dir = os.path.join(cfg['jobs_dir'], args.sample_id)
    os.makedirs(work_dir, exist_ok=True)
    dataset = SampleScanner(cfg).get_dataset(args.sample_id)
    sample_contamination_check = VerifyBamId(dataset, work_dir, args.bam_File)
    sample_contamination_check.start()
    exit_status = sample_contamination_check.join()
    return exit_status


def run_gender_validation(args):
    if args.work_dir:
        work_dir = args.work_dir
    else:
        work_dir = os.path.join(cfg['jobs_dir'], args.sample_id)
    os.makedirs(work_dir, exist_ok=True)
    dataset = SampleScanner(cfg).get_dataset(args.sample_id)
    s = GenderValidation(dataset, work_dir, args.vcf_file)
    s.start()
    return s.join()


if __name__ == '__main__':
    sys.exit(main())
