from .script_writer import ClusterWriter


class PBSWriter(ClusterWriter):
    """Writes a Bash script runnable on PBS"""
    suffix = '.pbs'
    array_index = 'PBS_ARRAY_INDEX'

    def _write_header(self, cpus, mem, job_name, queue, walltime=None, jobs=1):
        """Write a base PBS header. If multiple jobs, split them into a job array."""
        self.write_lines(
            '#!/bin/bash\n',
            '#PBS -l ncpus=%s,mem=%sgb' % (cpus, mem),
            '#PBS -q ' + queue,
            '#PBS -j ' + 'oe',
            '#PBS -o ' + self.log_file,
            '#PBS -W block=true'
        )
        if walltime:
            self.write_line('#PBS -l walltime=%s:00:00' % walltime)
        if job_name:
            self.write_line('#PBS -N ' + self._trim_field(job_name, 15))
        if jobs > 1:
            # specify a job array
            self.write_line('#PBS -J 1-' + str(jobs))
        self.write_line('cd $PBS_O_WORKDIR')
        self._line_break()
