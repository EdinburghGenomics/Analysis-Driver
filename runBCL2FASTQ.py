#!/opt/anaconda/bin/python

# create a PBS script to run BCL2FASTQ
fo = open("pbs/run.pbs", "wb")

fo.write("#!/bin/bash\n");
# walltime needed
fo.write( "#PBS -l walltime=00:10:00\n");

# PBS resources
fo.write("#PBS -l ncpus=4,mem=6gb\n");

# jobname
fo.write("#PBS -N bcl2fastq\n");

# queue name
fo.write("#PBS -q uv2000\n");

# input/output
fo.write("#PBS -j oe \n");

# output file name
fo.write("#PBS -o gccTest \n");

# working directory
fo.write("cd $PBS_O_WORKDIR\n");

# bash command to run bcl2fastq
fo.write("which gcc\n");

# close the PBS script
fo.close()

# Import the Python modules required for system operations
# To get environment variables
import os

# To run shell commands
import subprocess

# submit the job to the queue
p = subprocess.Popen(["qsub","pbs/run.pbs" ], stdout=subprocess.PIPE, stderr=subprocess.PIPE)

# gather output and errors
stdout, stderr = p.communicate()
print(stdout)

