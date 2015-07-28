__author__ = 'mwham'
import subprocess
import os
from .logger import AppLogger, NamedAppLogger
from . import fastq_handler
from analysis_driver import config


app_logger = NamedAppLogger('Util')


class AnalysisDriverError(Exception):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)


def localexecute(*args, stream=True, dry_run=False):
    app_logger.debug('Executing: ' + ' '.join(args))
    if dry_run:
        return 'dry_run', 'dry_run'

    proc = subprocess.Popen(list(args), stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    if stream:
        funclogger = NamedAppLogger(args[0])

        # log stdout while the process is running
        while proc.poll() is None:  # while no exit status
            line = proc.stdout.readline()
            if line:
                funclogger.info(line.decode('utf-8').rstrip('\n'))

        # log any remaining stdout after the process has finished
        for remaining_line in proc.stdout:
            funclogger.info(remaining_line.decode('utf-8').rstrip('\n'))

    else:
        out, err = proc.communicate()
        if type(out) is bytes:
            out = out.decode('utf-8')
        if type(err) is bytes:
            err = err.decode('utf-8')
        return out, err


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
    localexecute(
        'rsync',
        '-avu',
        '--exclude=Data',
        os.path.join(config.default['input_data_dir'], run_id),
        os.path.join(config.default['raw_dir'], run_id)
    )

