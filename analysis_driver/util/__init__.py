__author__ = 'mwham'
import subprocess
import logging
import os
from .logger import AppLogger


class AnalysisDriverError(Exception):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)


def localexecute(*args):
    proc = subprocess.Popen(list(args), stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    out, err = proc.communicate()
    return out.decode('utf-8'), err.decode('utf-8')


def shellexecute(*args):
    proc = subprocess.Popen(' '.join(args), stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
    out, err = proc.communicate()
    return out.decode('utf-8'), err.decode('utf-8')


def find_fastqs(path):
    fastqs = []
    for root, dirs, files in os.walk(path):
        for f in files:
            if f.endswith('.fastq.gz'):
                fastqs.append(os.path.join(root, f))
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

