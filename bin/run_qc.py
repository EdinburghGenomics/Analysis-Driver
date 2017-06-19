import os
import logging
import argparse
from sys import path
from egcg_core import executor, clarity, util
from egcg_core.clarity import get_user_sample_name
from egcg_core.app_logging import logging_default as log_cfg
from egcg_core.rest_communication import get_document
path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from analysis_driver import quality_control as qc
from analysis_driver.config import default as cfg, load_config
from analysis_driver.exceptions import AnalysisDriverError
from analysis_driver.dataset import NoCommunicationDataset, RunDataset
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
    geno_val.add_argument('--check_samples', nargs='*')
    geno_val.set_defaults(func=run_genotype_validation)

    sp_contamination = subparsers.add_parser('species_contamination_check')
    sp_contamination.add_argument('--fastq_files', required=True, nargs='+', help='fastq file pair')
    sp_contamination.set_defaults(func=run_species_contamination_check)

    sample_contamination = subparsers.add_parser('sample_contamination_check')
    sample_contamination.add_argument('--bam_file', required=True)
    sample_contamination.set_defaults(func=run_sample_contamination_check)

    gender_val = subparsers.add_parser('gender_validation')
    gender_val.add_argument('-v', '--vcf_file', dest='vcf_file', type=str, help='vcf file used to detect gender')
    gender_val.set_defaults(func=run_gender_validation)

    median_cov = subparsers.add_parser('median_coverage')
    median_cov.add_argument('--bam_file', required=True, help='the fastq file pairs')
    median_cov.set_defaults(func=median_coverage)

    contam_blast = subparsers.add_parser('contamination_blast')
    contam_blast.add_argument('--fastq_file', required=True, nargs='+', help='fastq file to check for contamination')
    contam_blast.set_defaults(func=contamination_blast)

    relatedness_parser = subparsers.add_parser('relatedness')
    relatedness_parser.add_argument('--gvcfs', required=True, nargs='+')
    relatedness_parser.add_argument('--reference', required=True)
    relatedness_parser.set_defaults(func=relatedness)

    bad_cycle_tile_parser = subparsers.add_parser('bad_cycle_tile')
    bad_cycle_tile_parser.add_argument('--window_size', type=int, default=50)
    bad_cycle_tile_parser.add_argument('--tile_quality_threshold', type=int, default=20)
    bad_cycle_tile_parser.add_argument('--cycle_quality_threshold', type=int, default=20)
    bad_cycle_tile_parser.set_defaults(func=detect_bad_cycle_tile_in_run)

    peddy_parser = subparsers.add_parser('peddy')
    peddy_parser.add_mutually_exclusive_group('--samples')
    peddy_parser.add_mutually_exclusive_group('--projects')
    peddy_parser.add_argument('--reference')
    peddy_parser.set_defaults(func=peddy)

    return parser.parse_args()


def run_genotype_validation(dataset, args):
    # Get the sample specific config
    cfg.merge(cfg['sample'])
    os.makedirs(os.path.join(cfg['jobs_dir'], dataset.name), exist_ok=True)

    sample_output_dir = os.path.join(cfg['output_dir'], args.project_id, dataset.name)
    genotype_vcfs = util.find_files(sample_output_dir, '*_genotype_validation.vcf.gz')
    if not genotype_vcfs:
        fastq_files = util.find_files(sample_output_dir, '*_R?.fastq.gz')
        genotype_vcf = None
    else:
        genotype_vcf = genotype_vcfs[0]
        fastq_files = []

    geno_val = qc.GenotypeValidation(
        dataset=dataset,
        fastq_files=sorted(fastq_files),
        vcf_file=genotype_vcf,
        check_neighbour=args.check_neighbour,
        check_project=args.check_project,
        list_samples=args.check_samples
    )
    geno_val.run()

    user_sample_id = get_user_sample_name(dataset.name)
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
    species_contamination_check = qc.ContaminationCheck(dataset=dataset, fastq_files=sorted(args.fastq_files))
    species_contamination_check.run()

    species_name = clarity.get_species_from_sample(dataset.name)
    fastqscreen_result = parse_fastqscreen_file(species_contamination_check.fastqscreen_expected_outfiles, species_name)
    print(fastqscreen_result)


def run_sample_contamination_check(dataset, args):
    os.makedirs(os.path.join(cfg['jobs_dir'], dataset.name), exist_ok=True)
    v = qc.VerifyBamID(dataset=dataset, bam_file=args.bam_file)
    v.run()


def run_gender_validation(dataset, args):
    os.makedirs(os.path.join(cfg['jobs_dir'], dataset.name), exist_ok=True)
    g = qc.GenderValidation(dataset=dataset, vcf_file=args.vcf_file)
    g.run()


def median_coverage(dataset, args):
    os.makedirs(os.path.join(cfg['jobs_dir'], dataset.name), exist_ok=True)
    s = qc.SamtoolsDepth(dataset=dataset, bam_file=args.bam_file)
    s.run()


def contamination_blast(dataset, args):
    os.makedirs(os.path.join(cfg['jobs_dir'], dataset.name), exist_ok=True)
    b = qc.ContaminationBlast(dataset=dataset, fastq_file=args.fastq_file)
    b.run()


def relatedness(dataset, args):
    os.makedirs(os.path.join(cfg['jobs_dir'], dataset.name), exist_ok=True)
    r = qc.Relatedness(dataset=dataset)
    r.run()
    o = qc.ParseRelatedness(dataset=dataset, parse_method='parse_relatedness')
    o.run()

def detect_bad_cycle_tile_in_run(dataset, args):
    cfg.merge(cfg['run'])
    dataset = RunDataset(args.dataset_name)
    d = qc.BadTileCycleDetector(
        dataset=dataset,
        window_size=args.window_size,
        tile_quality_threshold=args.tile_quality_threshold,
        cycle_quality_threshold=args.cycle_quality_threshold
    )
    bad_tiles = d.detect_bad_tile()
    bad_cycle = d.detect_bad_cycle()
    for lane in sorted(set(list(bad_tiles) + list(bad_cycle))):
        print('Lane %s' % lane)
        if lane in bad_cycle:
            print('Bad cycles are: ' + ', '.join([str(c) for c in bad_cycle[lane]]))
        if lane in bad_tiles:
            print('Bad tiles are: ' + ', '.join([str(c) for c in bad_tiles[lane]]))

def get_all_project_gvcfs(project_folder):
        gvcfs = []
        for path in os.walk(project_folder):
            gvcf = [i for i in path[2] if i.endswith('.g.vcf.gz')]
            if gvcf:
                gvcfs.append(os.path.join(path[0], ''.join(gvcf)))
        return gvcfs

def peddy(dataset, args):
    os.makedirs(os.path.join(cfg['jobs_dir'], dataset.name), exist_ok=True)
    all_gvcfs = []
    sample_ids = []
    if args.samples:
        cfg.merge(cfg['sample'])
        for sample_id in args.samples:
            sample_ids.append(sample_id)
            s = get_document(sample_id)
            project_id = s.get('project_id')
            gvcf_file = s.get('user_sample_id') + '.g.vcf.gz'
            gvcf = util.find_file(cfg['input_dir'], project_id, sample_id, gvcf_file)
            all_gvcfs.append(gvcf)
    elif args.projects:
        cfg.merge(cfg['project'])
        for project_id in args.projects:
            samples_for_project = clarity.get_sample_names_from_project(project_id)
            sample_ids.extend(samples_for_project)
            project_folder = util.find_file(cfg['input_dir'], project_id)
            all_gvcfs.extend(get_all_project_gvcfs(project_folder))


    g = qc.Genotype_gVCFs(dataset=dataset, GVCFs=all_gvcfs, reference=args.reference)
    g.run()
    p = qc.Peddy(dataset=dataset, ids=sample_ids)
    p.run()
    o = qc.ParseRelatedness(dataset=dataset, parse_method='parse_peddy')
    o.run()

if __name__ == '__main__':
    main()
