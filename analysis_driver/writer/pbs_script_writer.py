__author__ = 'mwham'
from .script_writer import ScriptWriter
from analysis_driver.config import default as cfg


class PBSWriter(ScriptWriter):
    """
    Writes a basic PBS submission script. Contains methods specific to PBS header statements and environment
    variables.
    """
    def __init__(self, script_name, walltime, cpus, mem, job_name, log_file, array=None, queue='uv2000'):
        """
        :param str script_name: As superclass
        :param int walltime: Desired walltime for the job
        :param int cpus: Number of cpus to allocate to the job
        :param int mem: Amount of memory to allocate to the job
        :param str job_name: A name to assign the job
        :param str log_file: A file to direct the job's stdout/stderr to
        :param int array: A number of jobs to submit in an array
        :param str queue: ID of the queue to submit to
        """
        super().__init__(script_name, array)
        self._write_header(walltime, cpus, mem, job_name, log_file, queue)
        self.info(
            'Written PBS header. Walltime %s, cpus %s, memory %s, job name %s, queue %s, array %s' % (
                walltime, cpus, mem, job_name, queue, array
            )
        )
        self.info('Job will write stdout to: ' + log_file)

    def start_array(self):
        self.write_line('case $PBS_ARRAY_INDEX in\n')

    def finish_array(self):
        self.write_line('*) echo "Unexpected PBS_ARRAY_INDEX: $PBS_ARRAY_INDEX"')
        self.write_line('esac')

    def _write_header(self, walltime, cpus, mem, job_name, log_file, queue):
        """
        Write a base PBS header using args from constructor.
        """
        self.write_line('#!/bin/bash\n')

        wt = self._write_pbs_param

        wt('-l', 'walltime=%s:00:00' % walltime)
        wt('-l', 'ncpus=%s,mem=%sgb' % (cpus, mem))
        if job_name:
            wt('-N', self._trim_field(job_name, 15))
        try:
            wt('-M', ','.join(cfg['notification_emails']))
            wt('-m', 'aeb')
        except KeyError:
            pass
        wt('-q', queue)
        wt('-j', 'oe')
        wt('-o', log_file)

        if self.array:
            wt('-J', ' 1-' + str(self.array) + '\n')

    def _write_pbs_param(self, flag, line):
        self.write_line('#PBS %s %s' % (flag, line))
