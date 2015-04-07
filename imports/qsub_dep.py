import distutils.spawn
import subprocess
import re

# PBS is installed in most of our computational resources
DEFAULT_QUEUING_SYSTEM = "PBS"

def qsub(args, queuing_system=DEFAULT_QUEUING_SYSTEM):
    """Submit job with ``qsub args`` and return the jobid.

    (REMEMBER: *args* is a list)
    """
    if queuing_system == "SLURM":
        base_cmd = "sbatch"
    else:
        base_cmd = "qsub"

    cmd = [base_cmd] + args
    print ">> " + " ".join(cmd)
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    output, errmsg = p.communicate()
    if p.returncode != 0:
        raise OSError(p.returncode, "command %r failed: %s" % (" ".join(cmd), errmsg))
    return get_jobid(output, queuing_system)


def get_jobid(s, queuing_system):
    """Get job id from the PBS
    """
    if queuing_system == "PBS":
        return s.strip()
    elif queuing_system == "GE":
        m =  re.search('Your job (?P<jobid>\d+) \("(?P<jobname>[^ "]+)"\)', s)
        if m:
            return m.group('jobid').strip()
    elif queuing_system == "SLURM":
        m =  re.search('Submitted batch job (?P<jobid>\d+)', s)
        if m:
            return m.group('jobid').strip()
    else:
        raise ValueError("Unknown queuing system %r" % queuing_system)
    return None

def dependent_job_args(jobid, queuing_system):
    templates = {'PBS':  ["-W", "depend=afterok:%s" % jobid],
                 'GE': ["-hold_jid", str(jobid)],
                 'SLURM': ["--dependency=afterok:%s" % jobid],
                 }
    return templates[queuing_system]

def qsub_dependents(args, jobid=None, queuing_system=DEFAULT_QUEUING_SYSTEM):
    """Submit jobs with *args*, possibly dependent on *jobid*.

    *args* is a list and contains everything on the ``qsub`` command
    line (although job dependency options should not be included).
    """
    if jobid is not None:
        new_args = dependent_job_args(jobid, queuing_system) + args
    else:
        new_args = args
    return qsub(new_args, queuing_system=queuing_system)
