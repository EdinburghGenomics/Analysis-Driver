import os
import logging
import argparse
import sys
from egcg_core import executor, util, rest_communication, archive_management
from egcg_core.app_logging import logging_default as log_cfg
from egcg_core import constants as c
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from analysis_driver.util import bash_commands
from analysis_driver.config import default as cfg, load_config
from analysis_driver.exceptions import PipelineError
from analysis_driver.dataset import RunDataset
from analysis_driver.report_generation import RunCrawler

app_logger = log_cfg.get_logger('Remove_phix')


def main():
    args = _parse_args()
    load_config()
    log_cfg.default_level = logging.DEBUG
    log_cfg.add_stdout_handler(logging.DEBUG)
    remove_phix(args.sample_id)


def _parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--sample_id', required=True)
    return parser.parse_args()


def remove_phix(sample_id):
    # Get the sample specific config
    cfg.merge(cfg['sample'])
    run_elements = rest_communication.get_documents('run_elements', where={'sample_id': sample_id})
    rundataset = RunDataset(run_elements[0][c.ELEMENT_RUN_NAME])

    rundataset.resolve_pipeline_and_toolset()
    job_dir = os.path.join(cfg['jobs_dir'], 'manual_' + sample_id)
    os.makedirs(job_dir, exist_ok=True)

    # Find all fastq files
    fastq_pairs = []
    for run_element in run_elements:
        run_dir = os.path.join(cfg['input_dir'], run_element[c.ELEMENT_RUN_NAME])
        fqs = util.find_fastqs(
            run_dir,
            run_element[c.ELEMENT_PROJECT_ID],
            run_element[c.ELEMENT_SAMPLE_INTERNAL_ID],
            run_element[c.ELEMENT_LANE]
        )
        fqs.sort()
        if run_element.get(c.ELEMENT_NB_READS_PASS_FILTER, 0) > 0:
            assert len(fqs) == 2
            fastq_pairs.append(fqs)

    # bwa alignment
    cmds = []
    for fqs in fastq_pairs:
        read_name_list = fqs[0][:-len('_R1_001.fastq.gz')] + '_phix_read_name.list'
        cmds.append(bash_commands.bwa_mem_phix(fqs[0], read_name_list))
    bwa_status = executor.execute(
        *cmds,
        job_name='phix_detection',
        working_dir=job_dir,
        cpus=16,
        mem=10,
        log_commands=False
    ).join()

    if bwa_status != 0:
        raise PipelineError('Bwa execution for sample %s failed with status %s' % (sample_id, bwa_status))

    # fastq filterer
    cmds = []
    for fqs in fastq_pairs:
        read_name_list = fqs[0][:-len('_R1_001.fastq.gz')] + '_phix_read_name.list'
        cmds.append(bash_commands.fastq_filterer(fqs, read_name_list))

    filterer_status = executor.execute(
        *cmds,
        prelim_cmds=[bash_commands.fq_filt_prelim_cmd()],
        job_name='fastq_filterer',
        working_dir=job_dir,
        cpus=18,
        mem=10
    ).join()
    if filterer_status != 0:
        raise PipelineError('Fastq_filterer execution for sample %s failed with status %s' % (sample_id, filterer_status))

    # fastqc
    fastqc_exec = executor.execute(
        *[bash_commands.fastqc(f) for f in [fq for fqs in fastq_pairs for fq in fqs]],
        job_name='fastqc',
        working_dir=job_dir,
        cpus=1,
        mem=2
    )

    # Seqtk fqchk
    seqtk_exec = executor.execute(
        *[bash_commands.seqtk_fqchk(f) for f in [fq for fqs in fastq_pairs for fq in fqs]],
        job_name='fqchk',
        working_dir=job_dir,
        cpus=1,
        mem=2,
        log_commands=False
    )

    # md5 calculation
    md5_exec = executor.execute(
        *[bash_commands.md5sum(f) for f in [fq for fqs in fastq_pairs for fq in fqs]],
        job_name='md5sum',
        working_dir=job_dir,
        cpus=1,
        mem=2,
        log_commands=False
    )

    fastqc_status = fastqc_exec.join()
    if fastqc_status != 0:
        raise PipelineError('Fastqc execution for sample %s failed with status %s' % (sample_id, fastqc_status))

    seqtk_status = seqtk_exec.join()
    if seqtk_status != 0:
        raise PipelineError('Seqtk execution for sample %s failed with status %s' % (sample_id, seqtk_status))

    md5_status = md5_exec.join()
    if md5_status != 0:
        raise PipelineError('Md5 execution for sample %s failed with status %s' % (sample_id, md5_status))

    for run_id in set(r[c.ELEMENT_RUN_NAME] for r in run_elements if r.get(c.ELEMENT_NB_READS_PASS_FILTER, 0) > 0):
        rd = RunDataset(run_id)
        run_dir = os.path.join(cfg['input_dir'], run_id)
        try:
            crawler = RunCrawler(rd, run_dir=run_dir, stage=RunCrawler.STAGE_MAPPING)
            crawler.send_data()
            # Ensure that the tape is up to date
            archive_management.archive_directory(run_dir)
        except Exception as e:
            app_logger.error(
                'Upload of metrics to reporting app or archiving for run %s encountered a %s exception: %s',
                run_id,
                e.__class__.__name__,
                e
            )


if __name__ == '__main__':
    sys.exit(main())
