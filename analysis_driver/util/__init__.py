__author__ = 'mwham'
import subprocess
import os
from .logger import AppLogger
from analysis_driver import config

class Util(AppLogger):
    pass


app_logger = Util()


class AnalysisDriverError(Exception):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)


def localexecute(*args, dry_run=False):
    app_logger.debug('Executing: ' + ' '.join(args))
    if not dry_run:
        proc = subprocess.Popen(list(args), stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        out, err = proc.communicate()
        # TODO: stream the stdout to the log instead
        return out.decode('utf-8'), err.decode('utf-8')
    else:
        return 'dry_run', 'dry_run'


def find_fastqs(path):
    fastqs = []
    for root, dirs, files in os.walk(path):
        for f in files:
            if f.endswith('.fastq.gz'):
                fastqs.append(os.path.join(root, f))
                # TODO: merge this into BCBioCSVWriter._find_fastqs()
    return fastqs


def setup_bcbio_run(bcbio, template, csv_file, run_dir, fastqs):
    out, err = localexecute(
        bcbio,
        '-w',
        'template',
        template,
        run_dir,
        csv_file,
        *fastqs
    )
    return out, err


def demultiplex_feedback(run_id):
    out, err = localexecute(
        'rsync',
        '-avu',
        '--exclude=Data',
        os.path.join(config.default['input_data_dir'], run_id),
        os.path.join(config.default['raw_dir'], run_id)
    )
    app_logger.info(out)
    if err:
        app_logger.warn(err)
