__author__ = 'mwham'
from .script_writer import ScriptWriter
from analysis_driver.config import default as cfg


class PBSWriter(ScriptWriter):
    """
    Writes a basic PBS submission script. Subclassed by BCL2FastqWriter, FastqcWriter and BCBioWriter.
    Initialises with self.script as an empty string, which is appended by self._write_line. This string is
    then saved to self.script_file by self.save.
    """
    def __init__(self, script_name, walltime, cpus, mem, job_name, log_file, array=None, queue='uv2000'):
        # TODO: pass dates, ints, etc. to constructor
        """
        :param str script_name: As superclass
        :param int walltime: Desired walltime for the job
        :param int cpus: Number of cpus to allocate to the job
        :param int mem: Amount of memory to allocate to the job
        :param str job_name: A name to assign the job
        :param str log_file: A file to direct the job's stdout/stderr to
        :param int array: A number of jobs to submit in an array
        :param queue: ID of the queue to submit to
        """
        super().__init__(script_name, array)
        self._write_header(walltime, cpus, mem, job_name, log_file, queue)
        self.info(
            'Written PBS header. Walltime: %s, cpus: %s, memory: %s, job name: %s, queue: %s' % (
                walltime, cpus, mem, job_name, queue
            )
        )
        self.info('Job will write stdout to: ' + log_file)

    def finish_array(self):
        self.write_line('*) echo "Unexpected PBS_ARRAY_INDEX: $PBS_ARRAY_INDEX"')
        self.write_line('esac')

    def _write_header(self, walltime, cpus, mem, job_name, log_file, queue):
        """
        Write a base PBS header using args from constructor.
        """
        wt = self.write_line
        wt('#!/bin/bash\n')

        wt('#PBS -l walltime=%s:00:00' % walltime)
        wt('#PBS -l ncpus=%s,mem=%sgb' % (cpus, mem))
        if job_name:
            wt('#PBS -N %s' % self._trim_field(job_name, 15))
        try:
            wt('#PBS -M ' + ','.join(cfg['notification_emails']))
            wt('#PBS -m aeb')
        except KeyError:
            pass
        wt('#PBS -q %s' % queue)  # queue name
        wt('#PBS -j oe')  # stdout/stderr
        wt('#PBS -o %s' % log_file)  # output file name
        if self.array:
            wt('#PBS -J 1-' + str(self.array))
        self.line_break()
