import os
import logging
import argparse
from sys import path
from egcg_core import executor, clarity, util
from egcg_core.clarity import get_user_sample_name
from egcg_core.app_logging import logging_default as log_cfg

path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from analysis_driver import quality_control as qc
from analysis_driver.config import default as cfg, load_config
from analysis_driver.exceptions import AnalysisDriverError
from analysis_driver.dataset import NoCommunicationDataset
from analysis_driver.util.bash_commands import rsync_from_to
from analysis_driver.reader.demultiplexing_parsers import parse_fastqscreen_file


def main():
    args = _parse_args()
    load_config()
    log_cfg.default_level = logging.DEBUG
    log_cfg.add_stdout_handler(logging.DEBUG)

    dataset = NoCommunicationDataset(args.dataset_name)
    os.makedirs(os.path.join(cfg['jobs_dir'], args.dataset_name), exist_ok=True)
    args.func(dataset, args)


def _parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--dataset_name', required=True)
    subparsers = parser.add_subparsers()

    geno_val = subparsers.add_parser('genotype_validation')
    geno_val.add_argument('--project_id', required=True)
    geno_val.add_argument('--check_neighbour', action='store_true', default=False)
    geno_val.add_argument('--check_project', action='store_true', default=False)
    geno_val.add_argument('--check_samples', nargs='*')
    geno_val.set_defaults(func=run_genotype_validation)

    sp_contamination = subparsers.add_parser('species_contamination_check')
    sp_contamination.add_argument('--fastq_files', required=True, nargs='+', help='fastq file pair')
    sp_contamination.set_defaults(func=run_species_contamination_check)

    sample_contamination = subparsers.add_parser('sample_contamination_check')
    sample_contamination.add_argument('--bam_file', required=True)
    sample_contamination.set_defaults(func=run_sample_contamination_check)

    gender_val = subparsers.add_parser('gender_validation')
    gender_val.add_argument('-v', '--vcf_file', dest='vcf_file', type=str, help='vcf file used to detect gender')
    gender_val.set_defaults(func=run_gender_validation)

    median_cov = subparsers.add_parser('median_coverage')
    median_cov.add_argument('--bam_file', required=True, help='the fastq file pairs')
    median_cov.set_defaults(func=median_coverage)

    contam_blast = subparsers.add_parser('contamination_blast')
    contam_blast.add_argument('--fastq_file', required=True, nargs='+', help='fastq file to check for contamination')
    contam_blast.set_defaults(func=contamination_blast)

    relatedness_parser = subparsers.add_parser('relatedness')
    relatedness_parser.add_argument('--gvcfs', required=True, nargs='+')
    relatedness_parser.add_argument('--reference', required=True)
    relatedness_parser.set_defaults(func=relatedness)

    return parser.parse_args()


def run_genotype_validation(dataset, args):
    # Get the sample specific config
    cfg.merge(cfg['sample'])

    sample_output_dir = os.path.join(cfg['output_dir'], args.project_id, dataset.name)
    genotype_vcfs = util.find_files(sample_output_dir, '*_genotype_validation.vcf.gz')
    if not genotype_vcfs:
        fastq_files = util.find_files(sample_output_dir, '*_R?.fastq.gz')
        genotype_vcf = None
    else:
        genotype_vcf = genotype_vcfs[0]
        fastq_files = []

    geno_val = qc.GenotypeValidation(
        dataset=dataset,
        fastq_files=sorted(fastq_files),
        vcf_file=genotype_vcf,
        check_neighbour=args.check_neighbour,
        check_project=args.check_project,
        list_samples=args.check_samples
    )
    geno_val.run()

    user_sample_id = get_user_sample_name(dataset.name)
    output_commands = []
    for f in [geno_val.seq_vcf_file, geno_val.validation_results.get(dataset.name)]:
        if f and os.path.exists(f):
            out_file = os.path.join(
                sample_output_dir,
                os.path.basename(f).replace(dataset.name, user_sample_id)
            )
            output_commands.append(rsync_from_to(f, out_file))

    exit_status = executor.execute(
        *output_commands,
        job_name='output_results',
        working_dir=geno_val.job_dir,
        cpus=1,
        mem=2
    ).join()

    if exit_status != 0:
        raise AnalysisDriverError('Copy of the results files to remote has failed')


def run_species_contamination_check(dataset, args):
    species_contamination_check = qc.ContaminationCheck(dataset=dataset, fastq_files=sorted(args.fastq_files))
    species_contamination_check.run()

    species_name = clarity.get_species_from_sample(dataset.name)
    fastqscreen_result = parse_fastqscreen_file(species_contamination_check.fastqscreen_expected_outfiles, species_name)
    print(fastqscreen_result)


def run_sample_contamination_check(dataset, args):
    v = qc.VerifyBamID(dataset=dataset, bam_file=args.bam_file)
    v.run()


def run_gender_validation(dataset, args):
    g = qc.GenderValidation(dataset=dataset, vcf_file=args.vcf_file)
    g.run()


def median_coverage(dataset, args):
    s = qc.SamtoolsDepth(dataset=dataset, bam_file=args.bam_file)
    s.run()


def contamination_blast(dataset, args):
    b = qc.ContaminationBlast(dataset=dataset, fastq_file=args.fastq_file)
    b.run()


def relatedness(dataset, args):
    r = qc.Relatedness(dataset=dataset, gvcf_files=args.gvcf_files,
                       reference=args.reference, project_id=dataset.name)
    r.run()


def peddy(args):
    pass


if __name__ == '__main__':
    main()
