__author__ = 'mwham'
import threading
import subprocess
import select  # asynchronous IO
import os.path
from analysis_driver.app_logging import AppLogger
from analysis_driver import writer
from analysis_driver.exceptions import AnalysisDriverError
from analysis_driver.config import default as cfg


class SimpleExecutor(AppLogger):
    def __init__(self, cmd):
        self.cmd = cmd
        self._validate_file_paths()
        self.proc = None

    def join(self):
        """
        Set self.proc to a Popen and start.
        :rtype: tuple[bytes, bytes]
        :raises: AnalysisDriverError on any exception
        """
        try:
            out, err = self._process().communicate()
            for line in out.decode('utf-8').split('\n'):
                self.info(line)
            for line in err.decode('utf-8').split('\n'):
                self.error(line)
            return self.proc.poll()
        except Exception as e:
            raise AnalysisDriverError('Command failed: ' + self.cmd) from e

    def start(self):
        raise NotImplementedError

    def _process(self):
        """
        Translate self.cmd to a subprocess. Override to manipulate how the process is run, e.g. with different
        resource managers.
        :rtype: subprocess.Popen
        """
        self.info('Executing: ' + self.cmd)
        self.proc = subprocess.Popen(self.cmd.split(' '), stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        return self.proc

    def _validate_file_paths(self):
        for arg in self.cmd.split(' '):
            if arg.startswith('/') and not os.path.exists(arg):
                self.debug('Could not find file: ' + arg + '. Will the executed command create it?')


class StreamExecutor(threading.Thread, SimpleExecutor):
    def __init__(self, cmd):
        """
        :param str cmd: A shell command to be executed
        """
        self.exception = None
        SimpleExecutor.__init__(self, cmd)
        threading.Thread.__init__(self)

    def join(self, timeout=None):
        """
        Ensure that both the thread and the subprocess have finished, and return self.proc's exit status.
        """
        super().join(timeout=timeout)
        if self.exception:
            self._stop()
            raise AnalysisDriverError('self.proc command failed: ' + self.cmd)
        return self.proc.wait()

    def run(self):
        try:
            self._stream_output()
        except Exception as e:
            self.exception = e

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


class ClusterExecutor(StreamExecutor):
    def __init__(self, cmds, **kwargs):
        """
        :param list cmds: Full path to a PBS script (for example)
        """
        prelim_cmds = kwargs.pop('prelim_cmds', None)
        w = writer.get_script_writer(jobs=len(cmds), **kwargs)
        w.write_jobs(cmds, prelim_cmds=prelim_cmds)
        super().__init__(w.script_name)

    def _process(self):
        """
        As the superclass, but with a qsub call to a PBS script.
        :rtype: subprocess.Popen
        """
        cmd = cfg.get('qsub', 'qsub') + ' ' + self.cmd

        self.info('Executing: ' + cmd)
        self.proc = subprocess.Popen(cmd.split(' '), stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        return self.proc


class ArrayExecutor(StreamExecutor):
    def __init__(self, cmds, stream):
        """
        :param list cmds:
        :param bool simple:
        """
        super().__init__(cmds)
        self.executors = []
        self.exit_statuses = []
        self.stream = stream
        for c in cmds:
            self.executors.append(StreamExecutor(c))

    def run(self):
        try:
            if self.stream:
                for e in self.executors:
                    e.start()
                for e in self.executors:
                    self.exit_statuses.append(e.join())
            else:
                for e in self.executors:
                    e.start()
                    self.exit_statuses.append(e.join())
        except Exception as err:
            self.exception = err

    def join(self, timeout=None):
        threading.Thread.join(self, timeout)
        if self.exception:
            self._stop()
            raise AnalysisDriverError('Commands failed: ' + str(self.exit_statuses))
        self.info('Exit statuses: ' + str(self.exit_statuses))
        return sum(self.exit_statuses)

    def _validate_file_paths(self):
        for cmd in self.cmd:
            for arg in cmd.split(' '):
                if arg.startswith('/') and not os.path.exists(arg):
                    self.debug('Could not find file: ' + arg + '. Will the executed command create it?')
