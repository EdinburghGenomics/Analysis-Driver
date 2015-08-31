__author__ = 'mwham'
import glob
import os.path
from analysis_driver import writer
from analysis_driver.app_logging import get_logger
from analysis_driver.executor import StreamExecutor
from analysis_driver.config import default as cfg

app_logger = get_logger('util')


def bcbio_prepare_samples(job_dir, sample_id, id_fastqs):
    # setup the BCBio merged csv file
    bcbio_csv_file = writer.write_bcbio_csv(job_dir, sample_id, id_fastqs)
    app_logger.info('Setting up BCBio samples from ' + bcbio_csv_file)

    merged_dir = os.path.join(job_dir, 'merged')
    return_code = _localexecute(
        cfg['bcbio_prepare_samples'],
        '--out',
        merged_dir,
        '--csv',
        bcbio_csv_file
    )
    # find the merged fastq files
    new_fastq_files = glob.glob(os.path.join(merged_dir, sample_id + '_R?.fastq.gz'))
    if return_code == 0:
        return new_fastq_files


def setup_bcbio_run(template, csv_file, run_dir, *fastqs):
    """
    Call localexecute to run 'bcbio -w template' on relevant input files.
    :param str bcbio: Path to the bcbio_nextgen.py executable
    :param str template: Path to the yaml file to use as the run template.
    :param str csv_file: Path to the csv sample file as generated by BCBioCSVWriter
    :param str run_dir: Path to the run folder
    :param str fastqs: Full paths to each input fastq file
    :return: None
    """
    app_logger.info('Setting up BCBio run')
    _localexecute(
        cfg['bcbio_nextgen'],
        '-w',
        'template',
        template,
        run_dir,
        csv_file,
        *fastqs
    )


def transfer_output_data(dataset):
    out_dir = os.path.join(cfg['output_dir'], dataset)
    if not os.path.isdir(out_dir):
        os.mkdir(out_dir)
    _localexecute(
        'rsync',
        '-avu',
        os.path.join(cfg['jobs_dir'], dataset),
        cfg['output_dir']
    )


def _localexecute(*args):
    executor = StreamExecutor(list(args))
    executor.start()
    return_code = executor.join()
    app_logger.info('Exit status: ' + str(return_code))
    return return_code
