import sys
import os
import argparse
import logging

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from analysis_driver.app_logging import logging_default as log_cfg
from analysis_driver import executor, util
from analysis_driver.exceptions import AnalysisDriverError
from analysis_driver.util.bash_commands import rsync_from_to, is_remote_path
from analysis_driver.config import default as cfg
from analysis_driver.constants import ELEMENT_NB_READS_CLEANED, ELEMENT_RUN_NAME, ELEMENT_PROJECT_ID, ELEMENT_LANE
from analysis_driver.dataset_scanner import SampleScanner
from analysis_driver.transfer_data import prepare_sample_data


log_cfg.set_log_level(logging.DEBUG)
log_cfg.add_handler(logging.StreamHandler(stream=sys.stdout))


def main():
    args = _parse_args()
    cfg.merge(cfg['sample'])
    fix_run_for_sample(args.sample_id)


def _parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--sample_id', required=True)

    return parser.parse_args()


def copy_sample_back(dataset, working_dir):
    """
    rsync the fastq files and qc to an the sample input_dir.
    :param Dataset dataset: A dataset object
    :param working_dir: The place where the fixed fastq files are
    """
    run_names = set()
    project_id = None
    for run_element in dataset.run_elements:
        if int(run_element.get(ELEMENT_NB_READS_CLEANED, 0)) > 0:
            run_names.add(run_element.get(ELEMENT_RUN_NAME))
            project_id = run_element.get(ELEMENT_PROJECT_ID)
    for run_name in run_names:
        fastq_dir = os.path.join(working_dir, run_name+'/', )
        copy_run_back(dataset.name, project_id,  run_name, fastq_dir)


def copy_run_back(sample_id, project_id, run_name, fastq_dir_to_transfer):
    dest_dir = os.path.join(cfg['input_dir'], run_name, 'fastq', project_id, sample_id)
    if not os.path.isdir(dest_dir) and not is_remote_path(dest_dir):
        raise AnalysisDriverError('Cannot find Destination directory for %s'%(sample_id))
    command = rsync_from_to(' '.join(fastq_dir_to_transfer), dest_dir)
    exit_status = executor.execute(
            [command],
            job_name='tf_bak',
            working_dir=os.path.join(cfg['jobs_dir'], sample_id)
        ).join()
    return exit_status


def fix_run_for_sample(sample_id):
    working_dir = os.path.join(cfg.query('jobs_dir'), sample_id)
    cfg.merge(cfg['sample'])
    scanner = SampleScanner(cfg)
    dataset = scanner.get_dataset(sample_id)
    fastq_files = prepare_sample_data(dataset)
    fastq_files.sort()
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

    copy_sample_back(dataset, working_dir)


if __name__ == '__main__':
    sys.exit(main())
