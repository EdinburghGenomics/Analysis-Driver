__author__ = 'mwham'
from .script_writer import ScriptWriter


class PBSWriter(ScriptWriter):
    """
    Writes a Bash script to be run on PBS. Contains methods specific to PBS header statements and environment
    variables.
    """
    suffix = '.pbs'

    def __init__(self, job_name, run_id, cpus, mem, walltime=None, jobs=1, **kwargs):
        """
        :param int jobs: Number of jobs to submit, in an array if needed
        """
        super().__init__(job_name, run_id, jobs, **kwargs)
        self._write_header(cpus, mem, job_name, self.queue, walltime, jobs)
        self.info(
            'Written PBS header. Walltime %s, cpus %s, memory %s, job name %s, queue %s, array %s' % (
                walltime, cpus, mem, job_name, self.queue, jobs
            )
        )

    def _start_array(self):
        self.write_line('case $PBS_ARRAY_INDEX in\n')

    def _finish_array(self):
        self.write_line('*) echo "Unexpected PBS_ARRAY_INDEX: $PBS_ARRAY_INDEX"')
        self.write_line('esac')

    def _write_header(self, cpus, mem, job_name, queue, walltime=None, jobs=1):
        """
        Write a base PBS header using args from constructor. If multiple jobs, split them into a job array.
        """
        self.write_line('#!/bin/bash\n')

        wt = self.write_line

        wt('#PBS -l ncpus=%s,mem=%sgb' % (cpus, mem))
        if walltime:
            wt('#PBS -l walltime=%s:00:00' % walltime)
        if job_name:
            wt('#PBS -N ' + self._trim_field(job_name, 15))
        wt('#PBS -q ' + queue)
        wt('#PBS -j ' + 'oe')
        wt('#PBS -o ' + self.log_file)
        wt('#PBS -W block=true')

        if jobs > 1:
            wt('#PBS -J 1-' + str(jobs))
        self._line_break()
