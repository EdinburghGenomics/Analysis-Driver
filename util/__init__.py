__author__ = 'mwham'
import subprocess


def localexecute(*args):
    proc = subprocess.Popen(list(args), stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    out, err = proc.communicate()

    print(out.decode('utf-8'))
    print(err.decode('utf-8'))

