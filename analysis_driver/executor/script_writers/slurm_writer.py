from . import ClusterWriter


class SlurmWriter(ClusterWriter):
    """Writes a Bash script runnable on Slurm"""
    suffix = '.slurm'
    array_index = 'SLURM_ARRAY_TASK_ID'

    def _write_header(self, cpus, mem, job_name, queue, walltime=None, jobs=1):
        """Write a base PBS header. If multiple jobs, split them into a job array."""
        self.write_lines(
            '#!/bin/bash\n',
            '#SBATCH --mem=%sg' % mem,
            '#SBATCH --cpus-per-task=%s' % cpus,
            '#SBATCH --partition=' + queue,
            '#SBATCH --output=' + self.log_file
        )
        if walltime:
            self.write_line('#SBATCH --time=%s:00:00' % walltime)
        if job_name:
            self.write_line('#SBATCH --job-name="%s"' % job_name)
        if jobs > 1:
            # specify a job array
            self.write_line('#SBATCH --array=1-' + str(jobs))
        self.write_line('cd ' + self.working_dir)
        self._line_break()
