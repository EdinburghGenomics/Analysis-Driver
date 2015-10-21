__author__ = 'mwham'
from glob import glob
import shutil
import os.path
import hashlib
from . import fastq_handler
from analysis_driver import writer, executor
from analysis_driver.app_logging import get_logger
from analysis_driver.config import default as cfg

app_logger = get_logger('util')


def bcbio_prepare_samples(job_dir, sample_id, fastqs, user_sample_id=None):
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
    cmd = (
        os.path.join(cfg['bcbio'], 'bin', 'bcbio_prepare_samples.py'),
        '--out',
        merged_dir,
        '--csv',
        bcbio_csv_file
    )
    return_code = executor.execute([' '.join(cmd)], env='local').join()
    merged_fastqs = glob(os.path.join(merged_dir, sample_id + '_R?.fastq.gz'))
    if merged_fastqs and return_code == 0:
        return merged_fastqs


def transfer_output_files(sample_id, output_dir, source_path_mapping):
    """
    :param str sample_id: a sample id, acting as the basename for the file
    :param str output_dir: where to send the output file
    :param dict source_path_mapping: a mapping describing where to find each type of output file
    :return: an exit status, non-zero if files weren not found
    """
    exit_status = 0
    for f in cfg['output_files']:
        app_logger.debug('Transferring ' + str(f))
        source_dir = source_path_mapping[f['type']]  # f['type'] = 'vcf'; source_dir = bcbio_source_dir
        base_name = f['name'].replace('*', sample_id)  # '*.vcf.gz' -> sample_id.vcf.gz
        source_file = os.path.join(source_dir, base_name)  # bcbio_source_dir/sample_id.vcf.gz
        rename_to = f.get('rename_to')  # '.g.vcf.gz'
        if rename_to:
            base_name = base_name.replace(f['name'][1:], rename_to)  # sample_id.vcf.gz -> sample_id.g.vcf.gz

        output_file = os.path.join(output_dir, base_name)  # output_dir/sample_id.g.vcf.gz
        app_logger.info('Looking for file: ' + source_file)

        if os.path.isfile(source_file):
            app_logger.info('Found file. Copying to ' + output_file + '.')
            shutil.copyfile(
                source_file,
                output_file
            )
            app_logger.debug('Generating md5 checksum')
            md5 = hashlib.md5()
            with open(output_file, 'rb') as g:
                chunk = g.read(8192)
                while chunk:
                    md5.update(chunk)
                    chunk = g.read(8192)

            with open(output_file + '.md5', 'w') as h:
                h.write(md5.hexdigest())
            app_logger.info('Done')

        else:
            app_logger.info('Expected output file not found.')
            exit_status += 1

    return exit_status
