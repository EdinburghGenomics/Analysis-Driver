import os
import csv
from os.path import join, isfile
from datetime import datetime

import illuminate
from bitstring import ReadError
from egcg_core import executor
from egcg_core.constants import ELEMENT_LANE, ELEMENT_SAMPLE_INTERNAL_ID, ELEMENT_LIBRARY_INTERNAL_ID, \
    ELEMENT_PROJECT_ID, ELEMENT_BARCODE
from egcg_core.util import str_join
from egcg_core.app_logging import AppLogger, logging_default as log_cfg
from analysis_driver.config import default as cfg
from . import bash_commands
app_logger = log_cfg.get_logger('util')


def bcbio_prepare_samples_cmd(job_dir, sample_id, fastqs, user_sample_id):
    """
    Call bcbio_prepare_samples with a csv sample file and a list of fastqs.
    :param str job_dir: Full path to the run folder
    :param str sample_id: Unique internal ID to assign to the samples
    :param list fastqs: Full paths to each input fastq file
    :param str user_sample_id: External sample ID for output filenames
    """
    # setup the BCBio merged csv file
    bcbio_csv_file = _write_bcbio_csv(job_dir, sample_id, fastqs, user_sample_id)
    app_logger.info('Setting up BCBio samples from ' + bcbio_csv_file)

    merged_dir = os.path.join(job_dir, 'merged')
    return str_join(
        os.path.join(cfg['tools']['bcbio'], 'bin', 'bcbio_prepare_samples.py'),
        '--out',
        merged_dir,
        '--csv',
        bcbio_csv_file,
        separator=' '
    )


def _write_bcbio_csv(run_dir, sample_id, fastqs, user_sample_id):
    """Write out a simple csv mapping fastq files to a sample id."""
    csv_file = os.path.join(run_dir, 'samples_' + sample_id + '.csv')
    app_logger.info('Writing BCBio sample csv ' + csv_file)

    with open(csv_file, 'w') as f:
        writer = csv.writer(f)
        writer.writerow(['samplename', 'description'])
        for fq in fastqs:
            writer.writerow([fq, user_sample_id])

    return csv_file


def generate_samplesheet(dataset, filename):
    all_lines = [
        '[Header]', 'Date, ' + datetime.now().strftime('%d/%m/%Y'), 'Workflow, Generate FASTQ Only', '',
        '[Settings]', 'Adapter, AGATCGGAAGAGCACACGTCTGAACTCCAGTCA',
        'AdapterRead2, AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT', '', '[Data]',
        'Lane,Sample_ID,Sample_Name,Sample_Project,index'
    ]
    for run_element in dataset.run_elements:
        all_lines.append(','.join([
            run_element[ELEMENT_LANE],
            run_element[ELEMENT_SAMPLE_INTERNAL_ID],
            run_element[ELEMENT_LIBRARY_INTERNAL_ID],
            run_element[ELEMENT_PROJECT_ID],
            run_element[ELEMENT_BARCODE]
        ]))
    with open(filename, 'w') as open_samplesheet:
        open_samplesheet.write('\n'.join(all_lines) + '\n')
