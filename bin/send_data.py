__author__ = 'mwham'
import sys
import os
import argparse

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from analysis_driver.reader import SampleSheet
from analysis_driver.report_generation.report_crawlers import SampleCrawler, RunCrawler


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--test', action='store_true')
    subparsers = p.add_subparsers()

    run_parser = subparsers.add_parser('run')
    run_parser.add_argument('run_id')
    run_parser.add_argument('--samplesheet')
    run_parser.add_argument('--conversion_stats', nargs='?', default=None)
    run_parser.set_defaults(func=run_crawler)

    sample_parser = subparsers.add_parser('sample')
    sample_parser.add_argument('project_id')
    sample_parser.add_argument('sample_id')
    sample_parser.add_argument('--input_dir')
    sample_parser.set_defaults(func=sample_crawler)

    args = p.parse_args()
    c = args.func(args)
    if args.test:
        print(c)
    else:
        c.send_data()


def run_crawler(args):
    for f in (args.samplesheet, args.conversion_stats):
        assert os.path.isfile(f), 'Missing file: ' + f
    return RunCrawler(args.run_id, SampleSheet(args.samplesheet), args.conversion_stats)


def sample_crawler(args):
    assert os.listdir(args.input_dir)
    return SampleCrawler(args.sample_id, args.project_id, args.input_dir)


if __name__ == '__main__':
    main()
