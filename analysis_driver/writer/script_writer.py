__author__ = 'mwham'
from analysis_driver.util import AppLogger


class ScriptWriter(AppLogger):
    """
    Writes a basic PBS submission script. Subclassed by BCL2FastqWriter, FastqcWriter and BCBioWriter.
    Initialises with self.script as an empty string, which is appended by self._write_line. This string is
    then saved to self.script_file by self.save.
    """
    def __init__(self, script_name, array=None):
        # TODO: pass dates, ints, etc. to constructor
        """
        :param str script_name: Desired full path to the pbs script to write
        :param int array: A number of jobs to submit in an array
        """
        self.info('Writing: ' + script_name)
        self.script_file = open(script_name, 'w')
        self.lines = []
        self.array = array
        if self.array:
            self.array_len = 0

    def write_line(self, line):
        self.lines.append(line)

    def write_array_cmd(self, job_number, cmd):
        self.write_line(str(job_number) + ') ' + cmd + '\n' + ';;')
        self.array_len += 1

    def finish_array(self):
        raise NotImplementedError

    def line_break(self):
        self.lines.append('')

    def save(self):
        """
        Save self.script to self.script_file. Also closes it. This is important!
        """
        for line in self.lines:
            self.script_file.write(line + '\n')
        self.script_file.close()
        self.info('Finished writing script')
        if self.array and self.array_len != self.array:
            self.error('Bad number of array jobs: %s written, %s expected' % (self.array_len, self.array))

    @staticmethod
    def _trim_field(field, max_length):
        if len(field) > max_length:
            return field[0:max_length]
        else:
            return field
