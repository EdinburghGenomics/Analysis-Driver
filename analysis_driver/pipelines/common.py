import os
import time
import shutil
import yaml
from egcg_core import executor, clarity, util
from analysis_driver.util import bash_commands, bcbio_prepare_samples_cmd
from egcg_core.app_logging import logging_default as log_cfg
from analysis_driver.config import output_files_config, default as cfg
from analysis_driver.report_generation.report_crawlers import SampleCrawler
from analysis_driver.transfer_data import output_sample_data, create_links_from_bcbio

app_logger = log_cfg.get_logger('common')


def link_results_files(sample_id, sample_dir, output_fileset):
    dir_with_linked_files = os.path.join(sample_dir, 'linked_output_files')
    os.makedirs(dir_with_linked_files, exist_ok=True)

    # Create the links from the bcbio output to one directory
    create_links_from_bcbio(
        sample_id,
        sample_dir,
        output_files_config.query(output_fileset),
        dir_with_linked_files
    )
    return dir_with_linked_files


def output_data(dataset, sample_dir, sample_id, dir_with_linked_files):
    exit_status = 0

    # upload the data to the rest API
    project_id = clarity.find_project_name_from_sample(sample_id)
    c = SampleCrawler(sample_id, project_id, dir_with_linked_files)
    c.send_data()

    # md5sum
    dataset.start_stage('md5sum')
    md5sum_exit_status = executor.execute(
        *[bash_commands.md5sum(os.path.join(dir_with_linked_files, f)) for f in os.listdir(dir_with_linked_files)],
        job_name='md5sum',
        working_dir=sample_dir,
        cpus=1,
        mem=2,
        log_commands=False
    ).join()
    dataset.end_stage('md5sum', md5sum_exit_status)

    exit_status += md5sum_exit_status

    # transfer output data
    dataset.start_stage('data_transfer')
    transfer_exit_status = output_sample_data(sample_id, dir_with_linked_files, cfg['output_dir'])
    dataset.end_stage('data_transfer', transfer_exit_status)
    exit_status += transfer_exit_status

    return exit_status


def bcbio_prepare_sample(job_dir, sample_id, fastq_files):
    """
    Merge the fastq files per sample using bcbio prepare sample
    """
    user_sample_id = clarity.get_user_sample_name(sample_id, lenient=True)
    cmd = bcbio_prepare_samples_cmd(job_dir, sample_id, fastq_files, user_sample_id=user_sample_id)

    exit_status = executor.execute(cmd, job_name='bcbio_prepare_samples', working_dir=job_dir,).join()
    sample_fastqs = util.find_files(job_dir, 'merged', user_sample_id + '_R?.fastq.gz')

    app_logger.info('bcbio_prepare_samples finished with exit status ' + str(exit_status))

    return sample_fastqs



def cleanup(dataset_name):
    exit_status = 0
    # wait for all the previous PBS steps to be done writing to the folder before cleaning it up
    time.sleep(120)
    job_dir = os.path.join(cfg['jobs_dir'], dataset_name)
    cleanup_targets = [job_dir]
    intermediates_dir = cfg.get('intermediate_dir')
    if intermediates_dir:
        cleanup_targets.append(os.path.join(intermediates_dir, dataset_name))

    for t in cleanup_targets:
        app_logger.info('Cleaning up ' + t)
        try:
            shutil.rmtree(t)
        except (OSError, FileNotFoundError, NotADirectoryError) as e:
            app_logger.error(str(e))
            app_logger.warning('Could not remove: ' + t)
            exit_status += 1

    return exit_status
