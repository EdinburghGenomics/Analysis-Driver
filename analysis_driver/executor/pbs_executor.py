__author__ = 'mwham'
import subprocess
import threading
from .executor import StreamExecutor


class PBSExecutor(StreamExecutor):
    def __init__(self, script, block=False):
        """
        :param str script: Full path to a PBS script
        :param bool block: Whether to run the job in blocking ('monitor') mode
        """
        self.cmd = script
        self.block = block
        threading.Thread.__init__(self)

    def _process(self):
        """
        Override to supply a qsub command with self.script
        :rtype: subprocess.Popen
        """
        cmd = ['qsub']

        if self.block:
            cmd.append('-W')
            cmd.append('block=true')

        cmd.append(self.cmd)

        return subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
