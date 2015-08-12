__author__ = 'mwham'
import threading
import subprocess
import select  # asynchronous IO
from analysis_driver.util.logger import AppLogger
from analysis_driver import driver


class Executor(AppLogger):
    def __init__(self, cmd):
        self.cmd = cmd

    def _process(self):
        """
        Translate self.cmd to a subprocess. Override to manipulate how the process is run, e.g. with different
        resource managers
        :rtype: subprocess.Popen
        """
        return subprocess.Popen(self.cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    def run(self):
        out, err = self._process().communicate()
        return out, err


class StreamExecutor(Executor, threading.Thread):
    def __init__(self, cmd):
        """
        :param list cmd: A field-separated command to be executed
        """
        self.cmd = cmd
        threading.Thread.__init__(self)

    def run(self):
        """
        Run self._process and stream its stdout and stderr
        :return:
        """
        proc = self._process()
        read_set = [proc.stdout, proc.stderr]

        while read_set:
            rlist, wlist, xlist = select.select(read_set, [], [])

            for stream, emit in ((proc.stdout, self.info), (proc.stderr, self.error)):
                if stream in rlist:
                    line = stream.readline().decode('utf-8').strip()
                    if line:
                        emit(line)
                    else:
                        stream.close()
                        read_set.remove(stream)


class ProcessTrigger(Executor, threading.Thread):
    def __init__(self, **kwargs):
        self.kwargs = kwargs
        threading.Thread.__init__(self)

    def run(self):
        driver.pipeline(**self.kwargs)
