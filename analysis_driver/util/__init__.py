import csv
import os
from datetime import datetime

import itertools
from egcg_core.app_logging import logging_default as log_cfg
from egcg_core.constants import ELEMENT_LANE, ELEMENT_SAMPLE_INTERNAL_ID, ELEMENT_LIBRARY_INTERNAL_ID, \
    ELEMENT_PROJECT_ID, ELEMENT_BARCODE
from egcg_core.exceptions import EGCGError
from egcg_core.util import str_join

from analysis_driver.config import default as cfg
from analysis_driver.reader.run_info import Reads
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


def find_all_fastq_pairs_for_lane(location, lane):
    fastqs = []
    for name, dirs, files in os.walk(location):
        fastqs.extend(os.path.join(name, f) for f in files if f.endswith('.fastq.gz') and '_L00%s_' % lane in f)
    if len(fastqs) % 2 != 0:
        raise EGCGError('Expected even number of fastq files in %s, found %s' % (location, len(fastqs)))
    fastqs.sort()
    app_logger.info('Found %s fastqs in %s for lane %s', len(fastqs), location, lane)
    return list(zip(*[iter(fastqs)] * 2))


def get_ranges(list_int):

    for a, b in itertools.groupby(enumerate(sorted(list_int)), lambda x: x[1] - x[0]):
        b = list(b)
        yield b[0][1], b[-1][1]


def convert_bad_cycle_in_trim(bad_cycle_list, run_info):
    if not bad_cycle_list:
        return None, None
    read1_length = Reads.num_cycles(run_info.reads.upstream_read)
    read2_length = Reads.num_cycles(run_info.reads.upstream_read)
    index_length = sum(run_info.reads.index_lengths)

    # split bad cycles between read1 and read2
    bad_cycle_list1 = [c for c in bad_cycle_list if c <= read1_length]
    bad_cycle_list2 = [c - (read1_length + index_length) for c in bad_cycle_list if c > read1_length + index_length]

    ranges1 = list(get_ranges(bad_cycle_list1))
    ranges2 = list(get_ranges(bad_cycle_list2))

    if ranges1 and ranges1[-1][1] == read1_length and ranges1[-1][0] < read1_length:
        trim_r1 = ranges1[-1][0]-1
    else:
        trim_r1 = None
    if ranges2 and ranges2[-1][1] == read2_length and ranges2[-1][0] < read2_length:
        trim_r2 = ranges2[-1][0]-1
    else:
        trim_r2 = None

    return trim_r1, trim_r2



