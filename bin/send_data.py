import os
import sys
import logging
import argparse
import json
from egcg_core.app_logging import logging_default as log_cfg

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from analysis_driver.reader import SampleSheet
from analysis_driver.report_generation.report_crawlers import SampleCrawler, RunCrawler
from analysis_driver.config import default as cfg
from analysis_driver.reader.run_info import RunInfo

log_cfg.default_level = logging.DEBUG
log_cfg.add_handler(logging.StreamHandler(stream=sys.stdout), logging.DEBUG)


def main():
    if 'run' not in sys.argv and 'sample' not in sys.argv:
        print("no mode specified - use either 'run' or 'sample'")
        return 1
    parent = argparse.ArgumentParser(add_help=False)
    parent.add_argument('--test', action='store_true')
    p = argparse.ArgumentParser()
    subparsers = p.add_subparsers()

    run_parser = subparsers.add_parser('run', parents = [parent])
    run_parser.add_argument('run_id')
    run_parser.add_argument('--samplesheet')
    run_parser.add_argument('--conversion_stats', nargs='?', default=None)
    run_parser.add_argument('--run_dir', help='e.g. jobs/<run_id>')
    run_parser.set_defaults(func=run_crawler)

    sample_parser = subparsers.add_parser('sample', parents = [parent])
    sample_parser.add_argument('project_id')
    sample_parser.add_argument('sample_id')
    sample_parser.add_argument('--input_dir')
    sample_parser.set_defaults(func=sample_crawler)

    args = p.parse_args()

    return args.func(args)


def run_crawler(args):
    cfg.merge(cfg['run'])
    if args.run_dir:
        run_dir = args.run_dir
    else:
        run_dir = os.path.join(cfg.query('output_dir'), args.run_id)
    if args.conversion_stats:
        conversion_stats = args.conversion_stats
    else:
        conversion_stats = os.path.join(run_dir, 'Stats', 'ConversionStats.xml')
    run_info = RunInfo(os.path.join(run_dir))
    if args.samplesheet:
        samplesheet = SampleSheet(args.samplesheet, has_barcode=run_info.mask.has_barcodes)
    else:
        samplesheet = SampleSheet(os.path.join(run_dir, 'SampleSheet_analysis_driver.csv'),
                                  has_barcode=run_info.mask.has_barcodes)

    c = RunCrawler(args.run_id, samplesheet, conversion_stats, run_dir)
    if args.test:
        print(
            json.dumps(
                {
                    'run_elements': list(c.barcodes_info.values()),
                    'unexpected_barcodes': list(c.unexpected_barcodes.values()),
                    'lanes': list(c.lanes.values()),
                    'runs': c.run,
                    'samples': list(c.libraries.values()),
                    'project': list(c.projects.values())
                },
                indent=4
            )
        )
    else:
        c.send_data()
    return 0


def sample_crawler(args):
    assert os.listdir(args.input_dir)
    c = SampleCrawler(args.sample_id, args.project_id, args.input_dir)
    if args.test:
        print(json.dumps({'samples': c.sample}, indent=4))
    else:
        c.send_data()
    return 0


if __name__ == '__main__':
    sys.exit(main())