__author__ = 'mwham'
from analysis_driver.util.logger import AppLogger


class PBSWriter(AppLogger):
    def __init__(self, pbs_name, walltime, cpus, mem, job_name, log_file, queue='uv2000'):
        self.log('Writing PBS file: ' + pbs_name)
        self.pbs_file = open(pbs_name, 'w')
        self.script = ''
        self._write_header(walltime, cpus, mem, job_name, log_file, queue)
        self.log(
            'Written PBS header. Walltime: %s, cpus: %s, memory: %s, job name: %s, queue: %s' % (
                walltime, cpus, mem, job_name, queue
            )
        )
        self.log('Job will write stdout to: ' + log_file)

    def write_line(self, line):
        self.script += line + '\n'

    def _write_header(self, walltime, cpus, mem, job_name, log_file, queue):
        self.write_line('#!/bin/bash\n')

        self.write_line('#PBS -l walltime=%s:00:00' % walltime)
        self.write_line('#PBS -l ncpus=%s,mem=%sgb' % (cpus, mem))
        if job_name:
            self.write_line('#PBS -N %s' % job_name)
        self.write_line('#PBS -q %s' % queue)  # queue name
        self.write_line('#PBS -j oe')  # input/output
        self.write_line('#PBS -o %s' % log_file)  # output file name
        self.write_line('\n')

    def save(self):
        self.pbs_file.write(self.script)
        print(self.pbs_file)
        self.pbs_file.close()
        self.log('Closed pbs file')
