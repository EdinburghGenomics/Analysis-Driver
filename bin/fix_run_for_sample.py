import glob
import sys
import os
import argparse
import logging
from ctypes import util

from analysis_driver.dataset_scanner import SampleScanner
from analysis_driver.transfer_data import prepare_sample_data

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from analysis_driver.quality_control.genotype_validation import GenotypeValidation
from analysis_driver.config import logging_default as log_cfg
from analysis_driver import executor
from analysis_driver.exceptions import AnalysisDriverError
from analysis_driver.util.bash_commands import rsync_from_to, is_remote_path
from analysis_driver.config import default as cfg
from analysis_driver.clarity import get_user_sample_name

log_cfg.default_level = logging.DEBUG
log_cfg.add_handler('stdout', logging.StreamHandler(stream=sys.stdout), logging.DEBUG)


def main():
    args = _parse_args()

    args.func(args)


def _parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--project_id', required=True)
    parser.add_argument('--sample_id', required=True)

    return parser.parse_args()


def fix_run_for_sample(sample_id):
    working_dir = os.path.join(cfg.query('jobs_dir'), sample_id)
    cfg.merge(cfg['sample'])
    scanner = SampleScanner(cfg)
    dataset = scanner.get_dataset(sample_id)
    fastq_files = prepare_sample_data(dataset)
    fastq_files.sorted()
    fastq_pairs = list(zip(*[iter(fastq_files)]*2))
    exit_status = executor.execute(
        [util.bash_commands.sickle_paired_end_in_place(fqs) for fqs in fastq_pairs],
        job_name='sickle_filter',
        working_dir=working_dir,
        cpus=1,
        mem=2,
        log_commands=False
    ).join()

    if exit_status:
        return exit_status

    # fastqc
    fastqc_executor = executor.execute(
        [util.bash_commands.fastqc(fq) for fq in fastq_files],
        job_name='fastqc',
        working_dir=working_dir,
        cpus=1,
        mem=2
    )

    # seqtk fqchk
    seqtk_fqchk_executor = executor.execute(
        [util.bash_commands.seqtk_fqchk(fq) for fq in fastq_files],
        job_name='fqchk',
        working_dir=working_dir,
        cpus=1,
        mem=2,
        log_commands=False
    )

    # md5sum
    md5sum_executor = executor.execute(
        [util.bash_commands.md5sum(fq) for fq in fastq_files],
        job_name='md5sum',
        working_dir=working_dir,
        cpus=1,
        mem=2,
        log_commands=False
    )
    fastqc_exit_status = fastqc_executor.join()
    seqtk_exit_status = seqtk_fqchk_executor.join()
    md5sum_exit_status = md5sum_executor.join()

    exit_status = exit_status + fastqc_exit_status + seqtk_exit_status + md5sum_exit_status

    if exit_status:
        return exit_status



if __name__ == '__main__':
    sys.exit(main())
