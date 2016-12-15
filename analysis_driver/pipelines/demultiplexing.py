import os
import time
import shutil
from egcg_core import executor, clarity, util
from analysis_driver import reader
from analysis_driver.pipelines.common import cleanup
from analysis_driver.util import bash_commands
from analysis_driver.dataset_scanner import RunDataset
from analysis_driver.exceptions import PipelineError, SequencingRunError
from egcg_core.app_logging import logging_default as log_cfg
from analysis_driver.config import  default as cfg
from analysis_driver.quality_control.lane_duplicates import WellDuplicates
from analysis_driver.reader.version_reader import write_versions_to_yaml
from analysis_driver.report_generation.report_crawlers import RunCrawler
from analysis_driver.transfer_data import prepare_run_data, output_run_data

app_logger = log_cfg.get_logger('demultiplexing')



def demultiplexing_pipeline(dataset):
    """
    :param RunDataset dataset: Dataset object
    :return: Exit status
    :rtype: int
    """
    exit_status = 0

    dataset.start_stage('transfer')
    input_run_folder = prepare_run_data(dataset)
    dataset.end_stage('transfer')

    run_id = os.path.basename(input_run_folder)
    job_dir = os.path.join(cfg['jobs_dir'], run_id)
    fastq_dir = os.path.join(job_dir, 'fastq')
    app_logger.info('Input run folder (bcl data source): ' + input_run_folder)
    app_logger.info('Fastq dir: ' + fastq_dir)
    app_logger.info('Job dir: ' + job_dir)

    run_info = reader.RunInfo(input_run_folder)
    dataset.start_stage('setup')
    reader.transform_sample_sheet(input_run_folder, remove_barcode=not run_info.mask.has_barcodes)
    sample_sheet = reader.SampleSheet(os.path.join(input_run_folder, 'SampleSheet_analysis_driver.csv'), has_barcode=run_info.mask.has_barcodes)
    if not sample_sheet.validate(run_info.mask):
        raise PipelineError('Validation failed. Check SampleSheet.csv and RunInfo.xml.')

    # Send the information about the run to the rest API
    crawler = RunCrawler(run_id, sample_sheet)
    crawler.send_data()
    dataset.end_stage('setup')

    # Need to sleep to make sure the LIMS has had enough time to update itself
    time.sleep(900)
    run_status = clarity.get_run(run_id).udf.get('Run Status')
    # TODO: catch bcl2fastq error logs instead
    if run_status != 'RunCompleted':
        app_logger.error('Run status is \'%s\'. Stopping.', run_status)
        raise SequencingRunError(run_status)

    # bcl2fastq
    mask = sample_sheet.generate_mask(run_info.mask)
    app_logger.info('bcl2fastq mask: ' + mask)  # e.g: mask = 'y150n,i6,y150n'

    dataset.start_stage('bcl2fastq')
    exit_status += executor.execute(
        bash_commands.bcl2fastq(input_run_folder, fastq_dir, sample_sheet.filename, mask),
        job_name='bcl2fastq',
        working_dir=job_dir,
        cpus=8,
        mem=32
    ).join()
    dataset.end_stage('bcl2fastq', exit_status)
    if exit_status:
        return exit_status

    # well duplicates
    # This could be executed at the same time as bcl2fastq but I need the fastq directory to exist
    well_dup_exec = WellDuplicates(dataset, job_dir, fastq_dir, input_run_folder)
    well_dup_exec.start()

    # Filter the adapter dimer from fastq with fastq_filtered
    dataset.start_stage('fastq_filterer')
    exit_status = executor.execute(
        *[bash_commands.fastq_filterer_and_pigz_in_place(fqs) for fqs in util.find_all_fastq_pairs(fastq_dir)],
        job_name='fastq_filterer',
        working_dir=job_dir,
        cpus=18,
        mem=10
    ).join()
    dataset.end_stage('fastq_filterer', exit_status)
    if exit_status:
        return exit_status

    exit_status += well_dup_exec.join()
    if exit_status:
        return exit_status

    # check file integrity
    dataset.start_stage('integrity_check')
    integrity_executor = executor.execute(
        *[bash_commands.gzip_test(fq) for fq in util.find_all_fastqs(fastq_dir)],
        job_name='integrity_check',
        working_dir=job_dir,
        cpus=1,
        mem=2
    )

    # fastqc
    dataset.start_stage('fastqc')
    fastqc_executor = executor.execute(
        *[bash_commands.fastqc(fq) for fq in util.find_all_fastqs(fastq_dir)],
        job_name='fastqc',
        working_dir=job_dir,
        cpus=1,
        mem=2
    )

    # seqtk fqchk
    dataset.start_stage('seqtk_fqchk')
    seqtk_fqchk_executor = executor.execute(
        *[bash_commands.seqtk_fqchk(fq) for fq in util.find_all_fastqs(fastq_dir)],
        job_name='fqchk',
        working_dir=job_dir,
        cpus=1,
        mem=2,
        log_commands=False
    )

    # md5sum
    dataset.start_stage('md5sum')
    md5sum_executor = executor.execute(
        *[bash_commands.md5sum(fq) for fq in util.find_all_fastqs(fastq_dir)],
        job_name='md5sum',
        working_dir=job_dir,
        cpus=1,
        mem=2,
        log_commands=False
    )

    integrity_exit_status = integrity_executor.join()
    dataset.end_stage('integrity_check', integrity_exit_status)
    if integrity_exit_status:
        return integrity_exit_status

    fastqc_exit_status = fastqc_executor.join()
    dataset.end_stage('fastqc', fastqc_exit_status)
    seqtk_fqchk_exit_status = seqtk_fqchk_executor.join()
    dataset.end_stage('seqtk_fqchk', seqtk_fqchk_exit_status)
    md5_exit_status = md5sum_executor.join()
    dataset.end_stage('md5sum', md5_exit_status)

    # Copy the Samplesheet Runinfo.xml run_parameters.xml to the fastq dir
    for f in ['SampleSheet.csv', 'SampleSheet_analysis_driver.csv', 'runParameters.xml',
              'RunInfo.xml', 'RTAConfiguration.xml']:
        shutil.copy2(os.path.join(input_run_folder, f), os.path.join(fastq_dir, f))
    if not os.path.exists(os.path.join(fastq_dir, 'InterOp')):
        shutil.copytree(os.path.join(input_run_folder, 'InterOp'), os.path.join(fastq_dir, 'InterOp'))

    # Find conversion xml file and adapter file, and send the results to the rest API
    conversion_xml = os.path.join(fastq_dir, 'Stats', 'ConversionStats.xml')
    adapter_trim_file = os.path.join(fastq_dir, 'Stats', 'AdapterTrimming.txt')

    if os.path.exists(conversion_xml) and os.path.exists(adapter_trim_file):
        app_logger.info('Found ConversionStats and AdaptorTrimming. Sending data.')
        crawler = RunCrawler(run_id, sample_sheet, adapter_trim_file=adapter_trim_file,
                             conversion_xml_file=conversion_xml, run_dir=fastq_dir)
        # TODO: review whether we need this
        json_file = os.path.join(fastq_dir, 'demultiplexing_results.json')
        crawler.write_json(json_file)
        crawler.send_data()
    else:
        app_logger.error('ConversionStats or AdaptorTrimming not found.')
        exit_status += 1

    write_versions_to_yaml(os.path.join(fastq_dir, 'program_versions.yaml'))
    dataset.start_stage('data_transfer')
    transfer_exit_status = output_run_data(fastq_dir, run_id)
    dataset.end_stage('data_transfer', transfer_exit_status)
    exit_status += transfer_exit_status + fastqc_exit_status + md5_exit_status

    if exit_status == 0:
        dataset.start_stage('cleanup')
        exit_status += cleanup(run_id)
        dataset.end_stage('cleanup', exit_status)
    return exit_status

