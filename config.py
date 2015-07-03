import os

home = os.getenv('HOME')

fastq = os.path.join(home, 'fastq')
jobs = os.path.join(home, 'jobs')

job_execution = None  # 'pbs' 'local' None

