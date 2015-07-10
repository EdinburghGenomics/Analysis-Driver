__author__ = 'mwham'
import subprocess


def localexecute(*args):
    print(*args)
    ###### Friday afternoon, get this *args parsing to work
    proc = subprocess.Popen(list(args), stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    out, err = proc.communicate()

    print(out)
    print(err)

