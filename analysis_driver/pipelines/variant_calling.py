import os
from egcg_core import executor, clarity
from analysis_driver.pipelines.common import link_results_files, output_data, cleanup
from analysis_driver.pipelines.qc_pipelines import _bam_file_production
from egcg_core.app_logging import logging_default as log_cfg
from analysis_driver.config import default as cfg
from analysis_driver.reader.version_reader import write_versions_to_yaml


app_logger = log_cfg.get_logger('variant_calling')


def _gatk_var_calling(dataset, species):

    sample_id = dataset.name
    gatk_run_dir = os.path.join(cfg['jobs_dir'], sample_id, 'gatk_var_calling')
    user_sample_id = clarity.get_user_sample_name(sample_id, lenient=True)
    os.makedirs(gatk_run_dir, exist_ok=True)

    def gatk_cmd(run_cls, input_bam, output, xmx=16, nct=16, ext=None):
        base_cmd = ('java -Xmx{xmx}m -XX:+UseSerialGC -Djava.io.tmpdir={tmpdir} -jar {gatk} -R {ref} '
                    '-I {input_bam} -T {run_cls} --read_filter BadCigar --read_filter NotPrimaryAlignment '
                    '-o {output} -l INFO -U LENIENT_VCF_PROCESSING ')

        if ext:
            base_cmd += ext
        if nct > 1:
            base_cmd += ' -nct %s' % nct

        return base_cmd.format(
            xmx=str(xmx * 1000),
            tmpdir=gatk_run_dir,
            gatk=cfg['tools']['gatk'],
            ref=cfg['references'][species]['fasta'],
            run_cls=run_cls,
            input_bam=input_bam,
            output=output
        )

    dbsnp = cfg['references'][species]['dbsnp']
    known_indels = cfg['references'][species].get('known_indels')

    basename = os.path.join(gatk_run_dir, user_sample_id)
    sorted_bam = os.path.join(cfg['jobs_dir'], sample_id, sample_id + '.bam')
    recal_bam = basename + '_recal.bam'
    output_grp = basename + '.grp'
    output_intervals = basename + '.intervals'
    indel_realigned_bam = basename + '_indel_realigned.bam'
    sample_gvcf = basename + '.g.vcf'
    sample_gvcfgz = sample_gvcf + '.gz'

    base_recal = executor.execute(
        gatk_cmd('BaseRecalibrator', sorted_bam, output_grp, xmx=48, ext='--knownSites ' + dbsnp),
        job_name='gatk_base_recal',
        working_dir=gatk_run_dir,
        cpus=16,
        mem=64
    ).join()

    print_reads = executor.execute(
        gatk_cmd('PrintReads', sorted_bam, recal_bam, xmx=48, ext=' -BQSR ' + output_grp),
        job_name='gatk_print_reads',
        working_dir=gatk_run_dir,
        cpus=16,
        mem=64
    ).join()

    realign_target_cmd = gatk_cmd('RealignerTargetCreator', recal_bam, output_intervals, nct=1)
    if known_indels:
        realign_target_cmd += ' --known ' + known_indels
    realign_target = executor.execute(
        realign_target_cmd,
        job_name='gatk_realign_target',
        working_dir=gatk_run_dir,
        mem=16
    ).join()

    realign_cmd = gatk_cmd(
        'IndelRealigner',
        recal_bam,
        indel_realigned_bam,
        nct=1,
        ext='-targetIntervals ' + output_intervals
    )
    if known_indels:
        realign_cmd += ' --knownAlleles ' + known_indels
    realign = executor.execute(
        realign_cmd,
        job_name='gatk_indel_realign',
        working_dir=gatk_run_dir,
        mem=16
    ).join()

    haplotype_cmd = gatk_cmd(
        'HaplotypeCaller',
        indel_realigned_bam,
        sample_gvcf,
        xmx=48,
        ext=('--pair_hmm_implementation VECTOR_LOGLESS_CACHING -ploidy 2 --emitRefConfidence GVCF '
             '--variant_index_type LINEAR --variant_index_parameter 128000 --dbsnp ' + dbsnp)
    )
    for annot in ('BaseQualityRankSumTest', 'FisherStrand', 'GCContent', 'HaplotypeScore', 'HomopolymerRun',
                  'MappingQualityRankSumTest', 'MappingQualityZero', 'QualByDepth', 'ReadPosRankSumTest',
                  'RMSMappingQuality', 'DepthPerAlleleBySample', 'Coverage', 'ClippingRankSumTest',
                  'DepthPerSampleHC'):
        haplotype_cmd += ' --annotation ' + annot
    haplotype = executor.execute(
        haplotype_cmd,
        job_name='gatk_haplotype_call',
        working_dir=gatk_run_dir,
        cpus=16,
        mem=64
    ).join()

    bgzip = executor.execute(
        '%s %s' % (cfg['tools']['bgzip'], sample_gvcf),
        job_name='bgzip',
        working_dir=gatk_run_dir,
        cpus=1,
        mem=8
    ).join()
    tabix = executor.execute(
        '%s -p vcf %s' % (cfg['tools']['tabix'], sample_gvcfgz),
        job_name='tabix',
        working_dir=gatk_run_dir,
        cpus=1,
        mem=8
    ).join()

    return base_recal + print_reads + realign_target + realign + haplotype + bgzip + tabix


def var_calling_pipeline(dataset, species):
    sample_id = dataset.name
    sample_dir = os.path.join(cfg['jobs_dir'], sample_id)

    exit_status = _bam_file_production(dataset, species)
    exit_status += _gatk_var_calling(dataset, species)

    # link the bcbio file into the final directory
    dir_with_linked_files = link_results_files(sample_id, sample_dir, 'gatk_var_calling')

    write_versions_to_yaml(os.path.join(dir_with_linked_files, 'program_versions.yaml'))

    exit_status += output_data(dataset, sample_dir, sample_id, dir_with_linked_files)

    if exit_status == 0:
        dataset.start_stage('cleanup')
        exit_status += cleanup(sample_id)
        dataset.end_stage('cleanup', exit_status)

    return exit_status

