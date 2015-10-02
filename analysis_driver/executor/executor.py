__author__ = 'mwham'
import threading
import subprocess
import select  # asynchronous IO
import os.path
from analysis_driver.app_logging import AppLogger
from analysis_driver.exceptions import AnalysisDriverError
from analysis_driver.config import default as cfg


class Executor(AppLogger):
    """
    Non-threaded execution
    """
    def __init__(self, cmd):
        self.cmd = cmd
        self._validate_file_paths()
        self.proc = None

    def run(self):
        """
        Set self.proc to a Popen and start.
        :rtype: tuple[bytes, bytes]
        :raises: AnalysisDriverError on any exception
        """
        try:
            return self._process().communicate()
        except Exception as e:
            raise AnalysisDriverError('Command failed: ' + self.cmd) from e

    def _process(self):
        """
        Translate self.cmd to a subprocess. Override to manipulate how the process is run, e.g. with different
        resource managers.
        :rtype: subprocess.Popen
        """
        self.info('Executing: ' + str(self.cmd))
        self.proc = subprocess.Popen(self.cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        return self.proc

    def _validate_file_paths(self):
        for arg in self.cmd:
            if arg.startswith('/') and not os.path.exists(arg):
                self.warn('Could not find file: ' + arg)


class StreamExecutor(Executor, threading.Thread):
    def __init__(self, cmd):
        """
        :param list cmd: A shell command to be executed
        """
        self.exception = None
        Executor.__init__(self, cmd)
        threading.Thread.__init__(self)

    def _stream_output(self):
        """
        Run self._process and log its stdout/stderr until an EOF.
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

    def join(self, timeout=None):
        """
        Ensure that both the thread and the subprocess have finished, and return self.proc's exit status.
        """
        super().join(timeout=timeout)
        if self.exception:
            raise self.exception
        try:
            return self.proc.wait()
        except Exception as e:
            raise AnalysisDriverError('self.proc command failed: ' + str(self.cmd)) from e

    def run(self):
        try:
            self._stream_output()
        except Exception as e:
            self.exception = e


class ClusterExecutor(StreamExecutor):
    def __init__(self, script, block=True):
        """
        :param str script: Full path to a PBS script (for example)
        :param bool block: Whether to run the job in blocking mode
        """
        super().__init__([script])
        self.block = block

    def _process(self):
        """
        As the superclass, but with a qsub call to a PBS script.
        :rtype: subprocess.Popen
        """
        if cfg['job_execution'] == 'pbs':
            cmd = self._pbs_cmd(self.block)
        else:
            cmd = ['sh']

        cmd.extend(self.cmd)
        self.info('Executing: ' + str(cmd))
        self.proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        return self.proc

    @staticmethod
    def _pbs_cmd(block):
        cmd = [cfg['qsub']]
        if block:
            cmd.append('-W')
            cmd.append('block=true')
        return cmd
