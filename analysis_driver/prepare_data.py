
__author__ = 'mwham'
import os
from time import sleep
from analysis_driver.dataset_scanner import DATASET_READY, DATASET_NEW
from analysis_driver import executor
from analysis_driver.app_logging import get_logger
from analysis_driver.config import default as cfg
from analysis_driver.report_generation.demultiplexing_report import ELEMENT_RUN_NAME, ELEMENT_LANE, ELEMENT_PROJECT
from analysis_driver.util.fastq_handler import find_fastqs
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
    rsync_cmd = 'rsync -aqu --size-only --partial %s %s' % (os.path.join(from_dir, dataset.name), to_dir)

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