__author__ = 'mwham'
import os.path
from analysis_driver.app_logging import AppLogger
from analysis_driver.config import default as cfg


class ScriptWriter(AppLogger):
    """
    Writes a basic job submission script. Subclassed by PBSWriter. Initialises with self.lines as an empty
    list, which is appended by self.write_line. This list is then saved line by line to self.script_file by
    self.save.
    """
    suffix = '.sh'

    def __init__(self, job_name, run_id, jobs=1):
        """
        :param str job_name: Desired full path to the pbs script to write
        :param int jobs: A number of jobs to submit in an array
        """
        self.script_name = os.path.join(cfg['jobs_dir'], run_id, job_name + self.suffix)
        self.log_file = os.path.join(cfg['jobs_dir'], run_id, job_name + '.log')
        self.queue = cfg['job_queue']
        self.info('Writing: ' + self.script_name)
        self.info('Log file: ' + self.log_file)
        self.lines = []
        self.job_total = jobs
        self.jobs_written = 0

    def write_line(self, line):
        self.lines.append(line)

    def write_array_cmd(self, job_number, cmd, log_file=None):
        """
        :param int job_number: The index of the job (i.e. which number the job has in the array)
        :param str cmd: The command to write
        """
        line = str(job_number) + ') ' + cmd
        if log_file:
            line += ' > ' + log_file
        line += '\n' + ';;'
        self.write_line(line)
        self.jobs_written += 1

    def start_array(self):
        self.write_line('case $JOB_INDEX in')

    def finish_array(self):
        self.write_line('*) echo "Unexpected JOBINDEX: $JOB_INDEX')
        self.write_line('esac')

    def line_break(self):
        self.lines.append('')

    def save(self):
        """
        Save self.lines to self.script_file. Also closes it. Always close it.
        """
        if self.job_total > 1 and self.job_total != self.jobs_written:
            self.critical(
                'Bad number of array jobs: %s written, %s expected' % (self.jobs_written, self.job_total),
                ValueError
            )
        script_file = open(self.script_name, 'w')

        for line in self.lines:
            script_file.write(line + '\n')
        script_file.close()
        self.info('Closed ' + self.script_name)

    @staticmethod
    def _trim_field(field, max_length):
        """
        Required for, e.g, name fields which break PBS if longer than 15 chars
        :param str field: The input string to trim
        :param int max_length: What to trim it to
        :return: field, trimmed to the specified length
        :rtype: str
        """
        if len(field) > max_length:
            return field[0:max_length]
        else:
            return field
