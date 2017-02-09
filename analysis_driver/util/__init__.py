import os
import csv
from os.path import join, isfile
import illuminate
from analysis_driver.config import default as cfg
from analysis_driver.exceptions import AnalysisDriverError
from egcg_core.app_logging import logging_default as log_cfg
from egcg_core.executor import ArrayExecutor
from egcg_core.util import str_join
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
    """
    Write out a simple csv mapping fastq files to a sample id
    """
    csv_file = os.path.join(run_dir, 'samples_' + sample_id + '.csv')
    app_logger.info('Writing BCBio sample csv ' + csv_file)

    with open(csv_file, 'w') as f:
        writer = csv.writer(f)
        writer.writerow(['samplename', 'description'])
        for fq in fastqs:
            writer.writerow([fq, user_sample_id])

    return csv_file


def check_bcl_files(run_dir, run_info, validation_log):
    tile_ids = run_info.tiles
    run_info_ncycles = sum(int(e.attrib['NumCycles']) for e in run_info.mask.reads)
    all_bcls = []

    validated_bcls = []
    if isfile(validation_log):
        with open(validation_log, 'r') as f:
            validated_bcls = [l.rstrip('\n') for l in f.readlines()]

    all_cycles = _all_cycles(run_dir)
    uniq_cycles = sorted(set(all_cycles), reverse=True)
    if all_cycles[-1] > run_info_ncycles:
        raise AnalysisDriverError(
            'Number of cycles (%s) disagrees with RunInfo (%s)' % (all_cycles[-1], run_info_ncycles)
        )

    last_completed_cycle = None
    for c in uniq_cycles:
        if all_cycles.count(c) == len(tile_ids):
            last_completed_cycle = c
            break
    if last_completed_cycle is None:
        raise AnalysisDriverError('Bad cycles')

    for c in range(last_completed_cycle):
        cycle_id = 'C' + str(c + 1) + '.1'
        for t in tile_ids:
            lane = t[2]  # s_1_1101.bcl.gz
            lane_id = 'L00' + lane
            all_bcls.append(join(lane_id, cycle_id, t + '.bcl.gz'))

    files_to_check = [f for f in all_bcls if f not in validated_bcls]

    basecalls_dir = os.path.join(run_dir, 'Data', 'Intensities', 'BaseCalls')
    app_logger.info('Checking %s files', len(files_to_check))
    e = ArrayExecutor(['gzip -t %s' % join(basecalls_dir, f) for f in files_to_check], stream=True)
    e.start()
    e.join()

    corrupted_files = []
    with open(validation_log, 'a') as f:
        # cmds and e.executors are in the same order, so zip produces corresponding bcls/exit statuses
        for bcl, exit_status in zip(files_to_check, e.exit_statuses):
            if exit_status == 0:
                f.write(bcl + '\n')
            else:
                corrupted_files.append(bcl)

    app_logger.info('Found %s corrupted files: %s', len(corrupted_files), corrupted_files)
    return corrupted_files


def _all_cycles(run_dir):
    return illuminate.InteropDataset(run_dir).ExtractionMetrics().data['cycle']
