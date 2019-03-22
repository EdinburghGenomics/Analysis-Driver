import os
import logging
import argparse
from sys import path
from egcg_core import executor, util, clarity
from egcg_core.app_logging import logging_default as log_cfg
from egcg_core.rest_communication import get_document
path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from analysis_driver import quality_control as qc
from analysis_driver.config import default as cfg, load_config
from analysis_driver.exceptions import AnalysisDriverError
from analysis_driver.dataset import NoCommunicationDataset, RunDataset, NoCommunicationSampleDataset
from analysis_driver.util.bash_commands import rsync_from_to
from analysis_driver.reader.demultiplexing_parsers import parse_fastqscreen_file


def main():
    args = _parse_args()
    load_config()
    log_cfg.default_level = logging.DEBUG
    log_cfg.add_stdout_handler(logging.DEBUG)

    dataset = NoCommunicationDataset(args.dataset_name)
    args.func(dataset, args)


def _parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--dataset_name', required=True)
    subparsers = parser.add_subparsers()

    geno_val = subparsers.add_parser('genotype_validation')
    geno_val.add_argument('--project_id', required=True)
    geno_val.add_argument('--check_neighbour', action='store_true', default=False)
    geno_val.add_argument('--check_project', action='store_true', default=False)
    geno_val.add_argument('--check_samples', nargs='*', default=[])
    geno_val.set_defaults(func=run_genotype_validation)

    sp_contamination = subparsers.add_parser('species_contamination_check')
    sp_contamination.add_argument('--fastq_files', required=True, nargs='+', help='fastq file pair')
    sp_contamination.set_defaults(func=run_species_contamination_check)

    sample_contamination = subparsers.add_parser('sample_contamination_check')
    sample_contamination.add_argument('--bam_file', required=True)
    sample_contamination.set_defaults(func=run_sample_contamination_check)

    sex_check = subparsers.add_parser('sex_check')
    sex_check.add_argument('-v', '--vcf_file', dest='vcf_file', type=str, help='vcf file used to detect sex')
    sex_check.set_defaults(func=run_sex_check)

    median_cov = subparsers.add_parser('median_coverage')
    median_cov.add_argument('--bam_file', required=True, help='the fastq file pairs')
    median_cov.set_defaults(func=median_coverage)

    contam_blast = subparsers.add_parser('contamination_blast')
    contam_blast.add_argument('--fastq_file', required=True, nargs='+', help='fastq file to check for contamination')
    contam_blast.set_defaults(func=contamination_blast)

    bad_cycle_tile_parser = subparsers.add_parser('bad_cycle_tile')
    bad_cycle_tile_parser.add_argument('--window_size', type=int, default=50)
    bad_cycle_tile_parser.add_argument('--tile_quality_threshold', type=int, default=20)
    bad_cycle_tile_parser.add_argument('--cycle_quality_threshold', type=int, default=20)
    bad_cycle_tile_parser.set_defaults(func=detect_bad_cycles_and_tiles)

    relatedness_parser = subparsers.add_parser('calculate_relatedness')
    data_type = relatedness_parser.add_mutually_exclusive_group()
    data_type.add_argument('--samples', nargs='+')
    data_type.add_argument('--projects', nargs='+')
    relatedness_parser.add_argument('--reference', required=True)
    relatedness_parser.add_argument('--method', required=True)
    relatedness_parser.set_defaults(func=calculate_relatedness)

    return parser.parse_args()


def run_genotype_validation(dataset, args):
    # Get the sample specific config
    dataset = NoCommunicationSampleDataset(dataset.name)
    dataset.resolve_pipeline_and_toolset()
    os.makedirs(os.path.join(cfg['jobs_dir'], dataset.name), exist_ok=True)

    sample_output_dir = os.path.join(cfg['sample']['output_dir'], args.project_id, dataset.name)
    genotype_vcfs = util.find_files(sample_output_dir, '*_genotype_validation.vcf.gz')
    if not genotype_vcfs:
        fq_pattern = os.path.join(sample_output_dir, '*_R?.fastq.gz')
        genotype_vcf = None
    else:
        genotype_vcf = genotype_vcfs[0]
        fq_pattern = []

    geno_val = qc.GenotypeValidation(
        dataset=dataset,
        fq_pattern=fq_pattern,
        vcf_file=genotype_vcf,
        check_neighbour=args.check_neighbour,
        check_project=args.check_project,
        list_samples=args.check_samples
    )
    geno_val.run()

    user_sample_id = clarity.get_user_sample_name(dataset.name)
    output_commands = []
    for f in [geno_val.seq_vcf_file, geno_val.validation_results.get(dataset.name)]:
        if f and os.path.exists(f):
            out_file = os.path.join(
                sample_output_dir,
                os.path.basename(f).replace(dataset.name, user_sample_id)
            )
            output_commands.append(rsync_from_to(f, out_file))

    exit_status = executor.execute(
        *output_commands,
        job_name='output_results',
        working_dir=geno_val.job_dir,
        cpus=1,
        mem=2
    ).join()

    if exit_status != 0:
        raise AnalysisDriverError('Copy of the results files to remote has failed')


def run_species_contamination_check(dataset, args):
    os.makedirs(os.path.join(cfg['jobs_dir'], dataset.name), exist_ok=True)
    f = qc.FastqScreen(dataset=dataset, fastq_files=sorted(args.fastq_files))
    f.run()

    fastqscreen_result = parse_fastqscreen_file(f.fastqscreen_expected_outfiles, dataset.species)
    print(fastqscreen_result)


def run_sample_contamination_check(dataset, args):
    os.makedirs(os.path.join(cfg['jobs_dir'], dataset.name), exist_ok=True)
    v = qc.VerifyBamID(dataset=dataset, bam_file=args.bam_file)
    v.run()


def run_sex_check(dataset, args):
    os.makedirs(os.path.join(cfg['jobs_dir'], dataset.name), exist_ok=True)
    g = qc.SexCheck(dataset=dataset, vcf_file=args.vcf_file)
    g.run()


def median_coverage(dataset, args):
    os.makedirs(os.path.join(cfg['jobs_dir'], dataset.name), exist_ok=True)
    s = qc.SamtoolsDepth(dataset=dataset, bam_file=args.bam_file)
    s.run()


def contamination_blast(dataset, args):
    os.makedirs(os.path.join(cfg['jobs_dir'], dataset.name), exist_ok=True)
    b = qc.Blast(dataset=dataset, fastq_file=args.fastq_file)
    b.run()


def detect_bad_cycles_and_tiles(dataset, args):
    dataset = RunDataset(args.dataset_name)
    d = qc.BadTileCycleDetector(
        dataset=dataset,
        window_size=args.window_size,
        tile_quality_threshold=args.tile_quality_threshold,
        cycle_quality_threshold=args.cycle_quality_threshold
    )
    bad_tiles = d.detect_bad_tiles()
    bad_cycle = d.detect_bad_cycles()
    for lane in sorted(set(list(bad_tiles) + list(bad_cycle))):
        print('Lane %s' % lane)
        if lane in bad_cycle:
            print('Bad cycles are: ' + ', '.join([str(c) for c in bad_cycle[lane]]))
        if lane in bad_tiles:
            print('Bad tiles are: ' + ', '.join([str(c) for c in bad_tiles[lane]]))


def get_all_project_gvcfs(project_folder):
        gvcfs = []
        for dirname, dirs, filenames in os.walk(project_folder):
            gvcfs.extend([os.path.join(dirname, f) for f in filenames if f.endswith('.g.vcf.gz')])
        return gvcfs


def calculate_relatedness(dataset, args):
    if not (args.samples or args.projects):
        raise AnalysisDriverError('Require either --samples or --projects parameter to be set')
    os.makedirs(os.path.join(cfg['jobs_dir'], dataset.name), exist_ok=True)
    all_gvcfs = []
    sample_ids = []
    if args.samples:
        for sample_id in args.samples:
            sample_ids.append(sample_id)
            s = get_document(sample_id)
            project_id = s.get('project_id')
            gvcf_file = s.get('user_sample_id') + '.g.vcf.gz'
            gvcf = util.find_file(cfg['sample']['input_dir'], project_id, sample_id, gvcf_file)
            all_gvcfs.append(gvcf)
    elif args.projects:
        for project_id in args.projects:
            samples_for_project = clarity.get_sample_names_from_project(project_id)
            sample_ids.extend(samples_for_project)
            project_folder = util.find_file(cfg['project']['input_dir'], project_id)
            all_gvcfs.extend(get_all_project_gvcfs(project_folder))

    g = qc.GenotypeGVCFs(dataset=dataset, gVCFs=all_gvcfs, reference=args.reference)
    g.run()
    if args.method not in ('peddy', 'relatedness'):
        raise AnalysisDriverError('Choose either "peddy" or "relatedness" as method')
    if args.method == 'peddy':
        p = qc.Peddy(dataset=dataset, ids=sample_ids)
        p.run()
        o = qc.ParseRelatedness(dataset=dataset, parse_method='parse_peddy', ids=sample_ids)
        o.run()
    elif args.method == 'relatedness':
        r = qc.Relatedness(dataset=dataset)
        r.run()
        o = qc.ParseRelatedness(dataset=dataset, parse_method='parse_relatedness', ids=sample_ids)
        o.run()

if __name__ == '__main__':
    main()
