__author__ = 'mwham'
from analysis_driver.util.logger import AppLogger


class PBSWriter(AppLogger):
    """
    Writes a basic PBS submission script. Subclassed by BCL2FastqWriter, FastqcWriter and BCBioWriter.
    Initialises with self.script as an empty string, which is appended by self._write_line. This string is
    then saved to self.script_file by self.save.
    """
    def __init__(self, pbs_name, walltime, cpus, mem, job_name, log_file, queue='uv2000'):
        # TODO: pass dates, ints, etc. to constructor
        """
        :param pbs_name: Desired full path to the pbs script to write
        :param walltime: Desired walltime for the job
        :param cpus: Number of cpus to allocate to the job
        :param mem: Amount of memory to allocate to the job
        :param job_name: A name to assign the job
        :param log_file: A file to direct the job's stdout/stderr to
        :param queue: ID of the queue to submit to
        """
        self.info('Writing PBS file: ' + pbs_name)
        self.pbs_file = open(pbs_name, 'w')
        self.script = ''
        self._write_header(walltime, cpus, mem, job_name, log_file, queue)
        self.info(
            'Written PBS header. Walltime: %s, cpus: %s, memory: %s, job name: %s, queue: %s' % (
                walltime, cpus, mem, job_name, queue
            )
        )
        self.info('Job will write stdout to: ' + log_file)

    def write_line(self, line):
        self.script += line + '\n'

    def _write_header(self, walltime, cpus, mem, job_name, log_file, queue):
        """
        Write a base PBS header using args from constructor.
        """
        self.write_line('#!/bin/bash\n')

        self.write_line('#PBS -l walltime=%s:00:00' % walltime)
        self.write_line('#PBS -l ncpus=%s,mem=%sgb' % (cpus, mem))
        if job_name:
            self.write_line('#PBS -N %s' % self._trim_field(job_name, 15))
        # TODO: -M <email_address>, -m aeb
        self.write_line('#PBS -q %s' % queue)  # queue name
        self.write_line('#PBS -j oe')  # stdout/stderr
        self.write_line('#PBS -o %s' % log_file)  # output file name
        self.write_line('\n')

    @staticmethod
    def _trim_field(field, max_length):
        if len(field) > max_length:
            return field[0:max_length]
        else:
            return field

    def save(self):
        """
        Save self.script to self.script_file. Also closes self.pbs_script. This is important!
        """
        self.pbs_file.write(self.script)
        self.pbs_file.close()
        self.info('Closed pbs file')


# TODO: less direct inheritance. subclasses should not be PBS-specific
# base PBS/SGE writers, etc
# helper bcl2fastq, etc. command generators
#
