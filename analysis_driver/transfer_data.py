import os
from egcg_core import clarity, util, archive_management
from analysis_driver.config import default as cfg
from egcg_core.app_logging import logging_default as log_cfg
from egcg_core.constants import ELEMENT_RUN_NAME, ELEMENT_LANE, ELEMENT_PROJECT_ID, ELEMENT_NB_READS_CLEANED

app_logger = log_cfg.get_logger(__name__)


def prepare_sample_data(dataset):
    app_logger.debug('Preparing dataset %s (%s)', dataset.name, dataset.dataset_status)
    fastqs = []
    for run_element in dataset.run_elements:
        if int(run_element.get(ELEMENT_NB_READS_CLEANED, 0)) > 0:
            fastqs.extend(_find_fastqs_for_sample(dataset.name, run_element))
    return fastqs


def _find_fastqs_for_sample(sample_id, run_element):
    local_fastq_dir = os.path.join(cfg['input_dir'], run_element.get(ELEMENT_RUN_NAME))
    app_logger.debug('Searching for fastqs in ' + local_fastq_dir)
    return util.find_fastqs(
        local_fastq_dir,
        run_element.get(ELEMENT_PROJECT_ID),
        sample_id,
        run_element.get(ELEMENT_LANE)
    )


def create_links_from_bcbio(sample_id, input_dir, output_cfg, link_dir):
    exit_status = 0
    user_sample_id = clarity.get_user_sample_name(sample_id, lenient=True)

    links = []

    for output_record in output_cfg.content.values():
        src_pattern = os.path.join(
            input_dir,
            os.path.join(*output_record.get('location', [''])),
            output_record['basename']
        ).format(sample_id=sample_id, user_sample_id=user_sample_id)

        sources = util.find_files(src_pattern)
        if sources:
            source = sources[-1]
            link_file = os.path.join(
                link_dir,
                output_record.get('new_name', os.path.basename(source))
            ).format(user_sample_id=user_sample_id)
            if not os.path.islink(link_file):
                os.symlink(source, link_file)
            links.append(link_file)
        else:
            app_logger.warning('No files found for pattern ' + src_pattern)
            if output_record.get('required', True):
                exit_status += 1
    if exit_status == 0:
        return links
    else:
        app_logger.error('link creation failed with exit status ' + str(exit_status))


def _output_data_and_archive(source_dir, output_dir,):
    exit_status = util.move_dir(source_dir, output_dir)
    if exit_status == 0:
        if archive_management.archive_directory(output_dir):
            return 0
        else:
            return 1
    return exit_status


def output_project_data(source_dir, project_id):
    return _output_data_and_archive(
        source_dir,
        os.path.join(cfg['project']['input_dir'], project_id, 'project_qc')
    )


def output_run_data(fastq_dir, run_id):
    """Retrieve and copy the fastq files to the output directory"""
    return _output_data_and_archive(
        fastq_dir,
        os.path.join(cfg['output_dir'], run_id)
    )


def output_sample_data(sample_id, source_dir, output_dir):
    project_id = clarity.find_project_name_from_sample(sample_id)
    output_dir = os.path.join(output_dir, project_id, sample_id)

    return _output_data_and_archive(
        source_dir.rstrip('/') + '/',
        output_dir.rstrip('/')
    )
