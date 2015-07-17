__author__ = 'mwham'
import subprocess
import logging
import os

def localexecute(*args):
    logging.info('Locally executing: ' + ' '.join(args))
    proc = subprocess.Popen(list(args), stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    out, err = proc.communicate()
    print(out.decode('utf-8'))
    print(err.decode('utf-8'))
    logging.info('Done')

def shellexecute(*args):
    logging.info('Locally executing through shell: ' + ' '.join(args))
    proc = subprocess.Popen(' '.join(args), stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
    out, err = proc.communicate()
    print(out.decode('utf-8'))
    print(err.decode('utf-8'))
    logging.info('Done')


def find_fastqs(path):
    fastqs = []
    for root, dirs, files in os.walk(path):
        for f in files:
            if f.endswith('.fastq.gz'):
                fastqs.append(os.path.join(root, f))
    return fastqs


class AnalysisDriverError(Exception):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)



