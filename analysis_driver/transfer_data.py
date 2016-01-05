__author__ = 'mwham'
import glob
import os
from time import sleep
from analysis_driver import executor, clarity
from analysis_driver.exceptions import AnalysisDriverError
from analysis_driver.writer.bash_commands import rsync_from_to, is_remote_path
from analysis_driver.app_logging import get_logger
from analysis_driver.config import default as cfg
from analysis_driver.report_generation.report_crawlers import ELEMENT_RUN_NAME, ELEMENT_LANE, ELEMENT_PROJECT
from analysis_driver.util import fastq_handler

app_logger = get_logger(__name__)


def prepare_run_data(dataset):
    """
    Decide whether to rsync a dataset to an intermediate dir and run driver.pipeline on it.
    :param Dataset dataset: A dataset object
    """
    app_logger.debug('Preparing dataset %s (%s)' % (dataset.name, dataset.dataset_status))
    if cfg.get('intermediate_dir'):
        _transfer_run_to_int_dir(
            dataset,
            cfg['input_dir'],
            cfg['intermediate_dir'],
            repeat_delay=int(cfg.get('tt_agent_delay', 120))
        )
        dataset_dir = cfg['intermediate_dir']
    else:
        dataset_dir = cfg['input_dir']

    return os.path.join(dataset_dir, dataset.name)


def prepare_sample_data(dataset):
    """
    Decide whether to rsync the fastq files to an intermediate dir just find them.
    :param Dataset dataset: A dataset object
    """
    app_logger.debug('Preparing dataset %s (%s)' % (dataset.name, dataset.dataset_status))
    fastqs = []

    for run_element in dataset.run_elements.values():
        fastqs.extend(_find_fastqs_for_sample(dataset.name, run_element))
    return fastqs


def _find_fastqs_for_sample(sample_id, run_element):
    run_id = run_element.get(ELEMENT_RUN_NAME)
    project_id = run_element.get(ELEMENT_PROJECT)
    lane = run_element.get(ELEMENT_LANE)

    local_fastq_dir = os.path.join(cfg['jobs_dir'], run_id, 'fastq')
    app_logger.debug('Searching for fastqs in ' + local_fastq_dir)
    fastqs = fastq_handler.find_fastqs(local_fastq_dir, project_id, sample_id, lane)
    if fastqs:
        return fastqs

    remote_fastq_dir = os.path.join(cfg['input_dir'], run_id, 'fastq')
    app_logger.debug('Searching for fastqs in ' + remote_fastq_dir)
    fastqs = fastq_handler.find_fastqs(remote_fastq_dir, project_id, sample_id, lane)
    if fastqs:
        return fastqs

    elif is_remote_path(remote_fastq_dir):
        pattern = os.path.join(remote_fastq_dir, project_id, sample_id, '*L00%s*.fastq.gz' % lane)

        # rsync the remote fastqs to a unique jobs dir
        rsync_cmd = rsync_from_to(pattern, os.path.join(cfg['jobs_dir'], sample_id, run_id))
        # TODO: try and parallelise this (although this avoids spamming the rdf server)
        exit_status = executor.execute([rsync_cmd], job_name='transfer_sample', run_id=sample_id).join()
        app_logger.info('Transfer complete with exit status ' + str(exit_status))

        app_logger.info('Searching again for fastqs in ' + local_fastq_dir)
        fastqs = glob.glob(os.path.join(cfg['jobs_dir'], sample_id, run_id, '*L00%s*.fastq.gz' % lane))

    if len(fastqs) != 2:
        raise AnalysisDriverError(
            '%s fastqs found for %s/%s/%s/L00%s' % (len(fastqs), run_id, project_id, sample_id, lane)
        )
    return fastqs


def _transfer_run_to_int_dir(dataset, from_dir, to_dir, repeat_delay, rsync_append_verify=True):
    exit_status = 0
    app_logger.info('Starting run transfer')

    rsync_cmd = rsync_from_to(
        os.path.join(from_dir, dataset.name),
        to_dir,
        append_verify=rsync_append_verify,
        exclude='Thumbnail_Images'
    )

    while not dataset.rta_complete():
        exit_status += executor.execute(
            [rsync_cmd],
            job_name='transfer_run',
            run_id=dataset.name
        ).join()
        sleep(repeat_delay)

    # one more rsync after the RTAComplete is created. After this, everything should be synced
    sleep(repeat_delay)
    exit_status += executor.execute([rsync_cmd], job_name='rsync', run_id=dataset.name).join()
    assert os.path.isfile(os.path.join(dataset.path, 'RTAComplete.txt'))
    app_logger.info('Transfer complete with exit status ' + str(exit_status))
    return exit_status


def create_links_from_bcbio(sample_id, intput_dir, output_config, link_dir, query_lims=True):
    exit_status = 0
    user_sample_id = None
    if query_lims:
        user_sample_id = clarity.get_user_sample_name(sample_id)
    if not user_sample_id:
        user_sample_id = sample_id
    os.makedirs(link_dir, exist_ok=True)

    links = []
    for output_record in output_config:
        src_pattern = os.path.join(
            intput_dir,
            os.path.join(*output_record['location']),
            output_record['basename']
        ).format(runfolder=sample_id, sample_id=user_sample_id)

        sources = glob.glob(src_pattern)
        if sources:
            source = sources[-1]
            link_file = os.path.join(
                link_dir,
                output_record.get('new_name', os.path.basename(source))
            ).format(sample_id=user_sample_id)
            if not os.path.islink(link_file):
                os.symlink(source, link_file)
            links.append(link_file)
        else:
            app_logger.warning('No files found for pattern ' + src_pattern)
            exit_status += 1
    if exit_status == 0:
        return links


def _output_data(source_dir, output_dir, run_id, rsync_append=True):
    if is_remote_path(output_dir):
        app_logger.info('output dir is remote')
        host, path = output_dir.split(':')
        ssh_cmd = 'ssh %s mkdir -p %s' % (host, path)
        exit_status = executor.execute([ssh_cmd], env='local', stream=False).join()
        if exit_status:
            raise AnalysisDriverError('Could not create remote output dir: ' + output_dir)

    else:
        os.makedirs(output_dir, exist_ok=True)

    command = rsync_from_to(source_dir, output_dir, append_verify=rsync_append)
    return executor.execute(
        [command],
        job_name='data_output',
        run_id=run_id
    ).join()


def output_run_data(fastq_dir, run_id):
    """Retrieve and copy the fastq files to the output directory"""
    return _output_data(fastq_dir, os.path.join(cfg['output_dir'], run_id), run_id)


def output_sample_data(sample_id, source_dir, output_dir, query_lims=True, rsync_append=True):
    if query_lims:
        project_id = clarity.find_project_from_sample(sample_id)
    else:
        project_id = 'proj_' + sample_id
    output_project = os.path.join(output_dir, project_id)
    output_sample = os.path.join(output_project, sample_id)

    return _output_data(
        source_dir.rstrip('/') + '/',
        output_sample.rstrip('/'),
        sample_id,
        rsync_append=rsync_append
    )