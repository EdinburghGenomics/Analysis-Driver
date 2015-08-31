__author__ = 'mwham'
from .script_writer import ScriptWriter
from analysis_driver.config import default as cfg


class PBSWriter(ScriptWriter):
    """
    Writes a basic PBS submission script. Contains methods specific to PBS header statements and environment
    variables.
    """
    suffix = '.pbs'

    def __init__(self, job_name, run_id, walltime, cpus, mem, jobs=1):
        """
        :param str job_name: A name to assign the job
        :param int walltime: Desired walltime for the job
        :param int cpus: Number of cpus to allocate to the job
        :param int mem: Amount of memory to allocate to the job
        :param int array: A number of jobs to submit in an array
        :param str queue: ID of the queue to submit to
        """
        super().__init__(job_name, run_id, jobs)
        self._write_header(walltime, cpus, mem, job_name, self.queue, jobs)
        self.info(
            'Written PBS header. Walltime %s, cpus %s, memory %s, job name %s, queue %s, array %s' % (
                walltime, cpus, mem, job_name, self.queue, jobs
            )
        )

    def start_array(self):
        self.write_line('case $PBS_ARRAY_INDEX in\n')

    def finish_array(self):
        self.write_line('*) echo "Unexpected PBS_ARRAY_INDEX: $PBS_ARRAY_INDEX"')
        self.write_line('esac')

    def _write_header(self, walltime, cpus, mem, job_name, queue, jobs):
        """
        Write a base PBS header using args from constructor.
        """
        self.write_line('#!/bin/bash\n')

        wt = self.write_line

        wt('#PBS -l walltime=%s:00:00' % walltime)
        wt('#PBS -l ncpus=%s,mem=%sgb' % (cpus, mem))
        if job_name:
            wt('#PBS -N ' + self._trim_field(job_name, 15))
        notification_emails = cfg.get('notification_emails')
        if notification_emails:
            wt('#PBS -M ' + ','.join(notification_emails))  # TODO: replace with notifications
            wt('#PBS -m aeb')
        wt('#PBS -q ' + queue)
        wt('#PBS -j ' + 'oe')
        wt('#PBS -o ' + self.log_file)

        if jobs > 1:
            wt('#PBS -J 1-' + str(jobs))
        self.line_break()
