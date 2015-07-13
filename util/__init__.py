__author__ = 'mwham'
import subprocess
import logging


def localexecute(*args):
    logging.debug('Executing: ' + str(args))
    proc = subprocess.Popen(list(args), stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    out, err = proc.communicate()

    print(out.decode('utf-8'))
    print(err.decode('utf-8'))


class AnalysisDriverError(Exception):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
