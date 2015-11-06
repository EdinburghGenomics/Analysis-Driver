__author__ = 'mwham'
import shutil
import os.path
import hashlib
from . import fastq_handler
from analysis_driver import writer, executor
from analysis_driver.app_logging import get_logger
from analysis_driver.config import default as cfg

app_logger = get_logger('util')


def bcbio_prepare_samples_cmd(job_dir, sample_id, fastqs, user_sample_id=None):
    """
    Call bcbio_prepare_samples with a csv sample file and a list of fastqs.
    :param str job_dir: Full path to the run folder
    :param str sample_id: An ID to assign to the samples
    :param list fastqs: Full paths to each input fastq file
    """
    # setup the BCBio merged csv file
    bcbio_csv_file = writer.write_bcbio_csv(job_dir, sample_id, fastqs, user_sample_id=user_sample_id)
    app_logger.info('Setting up BCBio samples from ' + bcbio_csv_file)

    merged_dir = os.path.join(job_dir, 'merged')
    return ' '.join(
        (
            os.path.join(cfg['bcbio'], 'bin', 'bcbio_prepare_samples.py'),
            '--out',
            merged_dir,
            '--csv',
            bcbio_csv_file
        )
    )


def transfer_output_file(source, dest):
    """
    :param str source:
    :param str dest:
    :return: exit status
    """
    app_logger.info('Transferring file: ' + source)
    if os.path.isfile(source):
        app_logger.debug('Found file. Transferring.')
        shutil.copyfile(source, dest)

        app_logger.debug('Generating md5 checksum')
        command = "md5sum %s > %s"%(source, source+'.md5')
        exit_status = executor.execute([command], env = 'local').join()
        md5_key = None
        with open(source + '.md5') as md5_f:
            md5_key, file_name = md5_f.readline().split()
        if md5_key:
            with open(dest + '.md5', 'w') as md5_f:
                md5_f.write('%s  %s'%(md5_key, os.path.basename(dest)))
        else:
            exit_status=1
        app_logger.debug('Done')
        return exit_status
    else:
        app_logger.warning('Could not find output file: %s'%source)
        return 1
