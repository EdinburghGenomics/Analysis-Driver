__author__ = 'mwham'
import subprocess
import logging
import os

from .logger import AppLogger

app_logger = logging.getLogger(__name__)

def localexecute(*args):
    app_logger.info('localexecute: ' + ' '.join(args))

    proc = subprocess.Popen(list(args), stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    out, err = proc.communicate()
    app_logger.info(out.decode('utf-8'))
    app_logger.info(err.decode('utf-8'))
    app_logger.info('Done')

def shellexecute(*args):
    app_logger.info('shellexecute: ' + ' '.join(args))

    proc = subprocess.Popen(' '.join(args), stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
    out, err = proc.communicate()
    app_logger.info(out.decode('utf-8'))
    app_logger.info(err.decode('utf-8'))
    app_logger.info('Done')

def find_fastqs(path):
    fastqs = []
    for root, dirs, files in os.walk(path):
        for f in files:
            if f.endswith('.fastq.gz'):
                fastqs.append(os.path.join(root, f))
    return fastqs

def setup_bcbio_run(bcbio, template, csv_file, run_dir, fastqs):
    localexecute(
        bcbio,
        '-w',
        'template',
        template,
        run_dir,
        csv_file,
        *fastqs
    )


class AnalysisDriverError(Exception):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
