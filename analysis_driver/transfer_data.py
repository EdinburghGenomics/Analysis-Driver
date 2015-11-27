import glob
from analysis_driver.writer.bash_commands import rsync_from_to

__author__ = 'mwham'
import os
from time import sleep
from analysis_driver.dataset_scanner import DATASET_READY, DATASET_NEW
from analysis_driver import executor, clarity
from analysis_driver.app_logging import get_logger
from analysis_driver.config import default as cfg
from analysis_driver.report_generation.run_report import ELEMENT_RUN_NAME, ELEMENT_LANE, ELEMENT_PROJECT
from analysis_driver.util.fastq_handler import find_fastqs
from analysis_driver.notification import default as ntf
from analysis_driver import writer

app_logger = get_logger('proctrigger')


def prepare_run_data(dataset):
    """
    Decide whether to rsync a dataset to an intermediate dir and run driver.pipeline on it.
    :param Dataset dataset: A dataset object
    """
    status = dataset.dataset_status
    exit_status = 0
    if cfg.get('intermediate_dir'):
        assert status in [DATASET_NEW, DATASET_READY], 'Invalid dataset status: ' + status
        exit_status = _transfer_to_int_dir(
            dataset,
            cfg['input_dir'],
            cfg['intermediate_dir'],
            repeat_delay=int(cfg.get('tt_agent_delay', 120))
        )
        dataset_dir = cfg['intermediate_dir']
    else:
        assert status in [DATASET_READY], 'Invalid dataset status: ' + status
        dataset_dir = cfg['input_dir']

    return os.path.join(dataset_dir, dataset.name)


def _transfer_to_int_dir(dataset, from_dir, to_dir, repeat_delay):
    """
    rsync -aqu --size-only --partial from_dir/dataset to_dir
    """
    exit_status = 0
    app_logger.info('Starting transfer')

    rsync_cmd = rsync_from_to(os.path.join(from_dir, dataset.name), to_dir)

    while dataset.dataset_status != DATASET_READY:
        exit_status += executor.execute([rsync_cmd], job_name='rsync', run_id=dataset.name, walltime=36).join()
        sleep(repeat_delay)

    # one more rsync after the RTAComplete is created. After this, everything should be synced
    sleep(repeat_delay)
    exit_status += executor.execute([rsync_cmd], job_name='rsync', run_id=dataset.name, walltime=36).join()
    assert os.path.isfile(os.path.join(dataset.path, 'RTAComplete.txt'))
    app_logger.info('Transfer complete with exit status ' + str(exit_status))
    return exit_status

def find_run_location(run_id):
    fastq_dir = os.path.join(cfg['jobs_dir'], run_id, 'fastq')
    if not os.path.isdir(fastq_dir):
        app_logger.debug(fastq_dir + 'does not exist')
        fastq_dir = os.path.join(cfg['output_dir'], 'runs', run_id)
    if not os.path.isdir(fastq_dir):
        app_logger.debug(fastq_dir + 'does not exist')
        return None
    return fastq_dir

def prepare_sample_data(dataset):
    """
    Decide whether to rsync the fastq files to an intermediate dir just find them.
    :param Dataset dataset: A dataset object
    """
    status = dataset.dataset_status
    exit_status = 0
    all_fastqs = []
    for run_element in dataset.run_elements.values():
        run_location = find_run_location(run_element.get(ELEMENT_RUN_NAME))
        all_fastqs.extend(find_fastqs(run_location, run_element.get(ELEMENT_PROJECT), dataset.name, run_element.get(ELEMENT_LANE)))

    return all_fastqs

def output_run_data(fastq_dir, run_id):
    """Retrieve and copy the fastq files to the output directory"""
    output_dir = cfg['output_dir']
    output_run_dir = os.path.join(output_dir, 'runs', run_id)
    command = rsync_from_to(fastq_dir, output_run_dir)
    return executor.execute([command], job_name='final_copy', run_id=run_id, walltime=36).join()


def output_sample_data(fastq_dir, run_id):
    """Retrieve and copy the fastq files to the output directory"""
    output_dir = cfg['output_dir']
    output_run_dir = os.path.join(output_dir, 'runs', run_id)
    command = rsync_from_to(fastq_dir, output_run_dir)
    return executor.execute([command], job_name='final_copy', run_id=run_id, walltime=36).join()


def output_sample_data(sample_id, intput_dir, output_dir, output_config):
    exit_status = 0
    
    project_id = clarity.find_project_from_sample(sample_id)
    user_sample_id = clarity.get_user_sample_name(sample_id)
    if not user_sample_id:
        user_sample_id = sample_id

    dir_with_linked_files = os.path.join(intput_dir,"linked_output_files")
    os.makedirs(dir_with_linked_files)

    list_of_link = []
    for output_record in output_config:
        src_pattern = os.path.join(
            intput_dir,
            os.path.join(*output_record['location']),
            output_record['basename']
        ).format(runfolder=sample_id, sample_id=user_sample_id)

        sources = glob(src_pattern)
        if sources:
            source = sources[-1]
            link_file = os.path.join(
                dir_with_linked_files,
                output_record.get('new_name', os.path.basename(source))
            ).format(sample_id=user_sample_id)
            os.symlink(source, link_file)
            list_of_link.append(link_file)
        else:
            app_logger.warning('No files found for pattern ' + src_pattern)
            exit_status += 1

    ntf.start_stage('md5sum')
    md5sum_exit_status = executor.execute(
        [writer.bash_commands.md5sum(f) for f in list_of_link],
        job_name='md5sum',
        run_id=sample_id,
        walltime=6,
        cpus=1,
        mem=2,
        log_command = False
    ).join
    ntf.end_stage('md5sum', md5sum_exit_status)

    exit_status += md5sum_exit_status
    output_loc = os.path.join(output_dir, project_id, sample_id)
    if not os.path.isdir(output_loc):
        os.makedirs(output_loc)

    command = rsync_from_to(dir_with_linked_files, output_loc)
    exit_status += executor.execute([command], job_name='final_copy', run_id=sample_id, walltime=36).join()

    #with open(os.path.join(output_loc, 'run_config.yaml'), 'w') as f:
    #    f.write(cfg.report())


    return exit_status
