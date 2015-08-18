__author__ = 'mwham'
import threading
import subprocess
import select  # asynchronous IO
from analysis_driver.util.logger import AppLogger
from analysis_driver.config import default as cfg


class Executor(AppLogger):
    def __init__(self, cmd):
        self.cmd = cmd
        self.proc = None

    def run(self):
        out, err = self._process().communicate()
        return out, err

    def _process(self):
        """
        Translate self.cmd to a subprocess. Override to manipulate how the process is run, e.g. with different
        resource managers
        :rtype: subprocess.Popen
        """
        self.proc = subprocess.Popen(self.cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        return self.proc


class StreamExecutor(Executor, threading.Thread):
    def __init__(self, cmd):
        """
        :param list cmd: A field-separated command to be executed
        """
        Executor.__init__(self, cmd)
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

    def join(self, timeout=None):
        super().join(timeout=timeout)
        return self.proc.poll()


class ClusterExecutor(StreamExecutor):
    def __init__(self, script, block):
        """
        :param str script: Full path to a PBS script
        :param bool block: Whether to run the job in blocking ('monitor') mode
        """
        super().__init__([script])
        self.block = block

    def _process(self):
        """
        Override to supply a qsub command
        :rtype: subprocess.Popen
        """
        cmd = []
        if cfg['job_execution'] == 'pbs':
            cmd = self._pbs_cmd(self.block)

        cmd.extend(self.cmd)
        self.proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        return self.proc

    @staticmethod
    def _pbs_cmd(block):
        cmd = ['qsub']
        if block:
            cmd.append('-W')
            cmd.append('block=true')
        return cmd
