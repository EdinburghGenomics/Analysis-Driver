import threading
import subprocess
import select  # asynchronous IO
import os.path
import shlex
from time import sleep
from analysis_driver.app_logging import AppLogger
from analysis_driver.exceptions import AnalysisDriverError
from analysis_driver.config import default as cfg
from . import script_writers


class Executor(AppLogger):
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
            for stream, emit in ((out, self.info), (err, self.error)):
                for line in stream.decode('utf-8').split('\n'):
                    emit(line)
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
        # TODO: explore how commands that contains bash specific construct can be run ie: command <(sub command)
        # will need to add shell=True which make some of the test fail
        self.proc = subprocess.Popen(shlex.split(self.cmd), stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        return self.proc

    def _validate_file_paths(self):
        for arg in self.cmd.split(' '):
            if arg.startswith('/') and not os.path.exists(arg):
                self.debug('Could not find file: ' + arg + '. Will the executed command create it?')


class StreamExecutor(threading.Thread, Executor):
    def __init__(self, cmd):
        """
        :param str cmd: A shell command to be executed
        """
        self.exception = None
        Executor.__init__(self, cmd)
        threading.Thread.__init__(self)

    def join(self, timeout=None):
        """
        Ensure that both the thread and the subprocess have finished, and return self.proc's exit status.
        """
        super().join(timeout=timeout)
        if self.exception:
            self._stop()
            self.error(self.exception.__class__.__name__ + ': ' + str(self.exception))
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


class ArrayExecutor(StreamExecutor):
    def __init__(self, cmds, stream):
        """
        :param cmds:
        :param bool stream: Whether to run all commands in parallel or one after another
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
            self.error(self.exception.__class__.__name__ + ': ' + str(self.exception))
            raise AnalysisDriverError('Commands failed: ' + str(self.exit_statuses))
        self.info('Exit statuses: ' + str(self.exit_statuses))
        return sum(self.exit_statuses)

    def _validate_file_paths(self):
        pass


class ClusterExecutor(AppLogger):
    script_writer = None
    finished_statuses = None
    unfinished_statuses = None

    def __init__(self, *cmds, prelim_cmds=None, **cluster_config):
        """
        :param list cmds: Full path to a job submission script
        """
        self.job_queue = cfg['job_queue']
        self.job_id = None
        w = self._get_writer(jobs=len(cmds), **cluster_config)
        if cfg.get('pre_job_source'):
            if not prelim_cmds:
                prelim_cmds = []
            prelim_cmds.append('source ' + cfg.get('pre_job_source'))
        w.write_jobs(cmds, prelim_cmds)
        qsub = cfg.query('tools', 'qsub', ret_default='qsub')
        self.cmd = qsub + ' ' + w.script_name

    def start(self):
        self.job_id = self._submit_job()
        self.info('Submitted "%s" as job %s' % (self.cmd, self.job_id))

    def join(self):
        sleep(10)
        while not self._job_finished():
            sleep(30)
        return self._job_exit_code()

    @classmethod
    def _get_writer(cls, job_name, working_dir, walltime=None, cpus=1, mem=2, jobs=1, log_commands=True):
        return cls.script_writer(job_name, working_dir, cfg['job_queue'], cpus, mem, walltime, jobs, log_commands)

    def _job_status(self):
        raise NotImplementedError

    def _job_exit_code(self):
        raise NotImplementedError

    def _submit_job(self):
        p = self._get_stdout(self.cmd)
        if p is None:
            raise AnalysisDriverError('Job submissions failed')
        return p

    def _job_finished(self):
        status = self._job_status()
        if status in self.unfinished_statuses:
            return False
        elif status in self.finished_statuses:
            return True
        self.debug('Bad job status: %s', status)

    def _get_stdout(self, cmd):
        p = subprocess.Popen(cmd.split(' '), stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        exit_status = p.wait()
        o, e = p.stdout.read(), p.stderr.read()
        self.debug('%s -> (%s, %s, %s)', cmd, exit_status, o, e)
        if exit_status:
            return None
        else:
            return o.decode('utf-8').strip()


class PBSExecutor(ClusterExecutor):
    unfinished_statuses = 'BEHQRSTUW'
    finished_statuses = 'FXM'
    script_writer = script_writers.PBSWriter

    def _qstat(self):
        h1, h2, data = self._get_stdout('qstat -x {j}'.format(j=self.job_id)).split('\n')
        return data.split()

    def _job_status(self):
        job_id, job_name, user, time, status, queue = self._qstat()
        return status

    def _job_exit_code(self):
        return self.finished_statuses.index(self._job_status())


class SlurmExecutor(ClusterExecutor):
    unfinished_statuses = ('RUNNING', 'RESIZING', 'SUSPENDED', 'PENDING')
    finished_statuses = ('COMPLETED', 'CANCELLED', 'FAILED', 'TIMEOUT', 'NODE_FAIL')
    script_writer = script_writers.SlurmWriter

    def _submit_job(self):
        # sbatch stdout: "Submitted batch job {job_id}"
        return super()._submit_job().split()[-1].strip()

    def _sacct(self, output_format):
        return self._get_stdout('sacct -n -j {j} -o {o}'.format(j=self.job_id, o=output_format))

    def _squeue(self):
        s = self._get_stdout('squeue -j {j} -o %T'.format(j=self.job_id))
        if not s or len(s.split('\n')) < 2:
            return None
        return sorted(set(s.split('\n')[1:]))

    def _job_status(self):
        state = self._squeue()
        if not state:
            state = self._sacct('State')
        return state

    def _job_exit_code(self):
        state, exit_code = self._sacct('State,ExitCode').split()
        if state == 'CANCELLED':  # cancelled jobs can still be exit status 0
            return 9
        return int(exit_code.split(':')[0])

