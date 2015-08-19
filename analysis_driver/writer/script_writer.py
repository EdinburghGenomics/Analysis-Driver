__author__ = 'mwham'
from analysis_driver.util import AppLogger


class ScriptWriter(AppLogger):
    """
    Writes a basic job submission script. Subclassed by PBSWriter. Initialises with self.lines as an empty
    list, which is appended by self.write_line. This list is then saved line by line to self.script_file by
    self.save.
    """
    def __init__(self, script_name, array=None):
        """
        :param str script_name: Desired full path to the pbs script to write
        :param int array: A number of jobs to submit in an array
        """
        self.script_name = script_name
        self.info('Writing: ' + self.script_name)
        self.lines = []
        self.array = array
        if self.array:
            self.array_len = 0

    def write_line(self, line):
        self.lines.append(line)

    def write_array_cmd(self, job_number, cmd):
        """
        :param int job_number: The index of the job (i.e. which number the job has in the array)
        :param str cmd: The command to write
        """
        self.write_line(str(job_number) + ') ' + cmd + '\n' + ';;')
        self.array_len += 1

    def start_array(self):
        raise NotImplementedError

    def finish_array(self):
        raise NotImplementedError

    def line_break(self):
        self.lines.append('')

    def save(self):
        """
        Save self.lines to self.script_file. Also closes it. Always close it.
        """
        if self.array and self.array_len != self.array:
            self.critical(
                'Bad number of array jobs: %s written, %s expected' % (self.array_len, self.array),
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
