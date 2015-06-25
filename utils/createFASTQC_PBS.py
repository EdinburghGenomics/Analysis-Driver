#!/opt/anaconda/bin/python

# Creates the PBS job file to run fastqc on all the fastq files in
# the specified directory or its subdirectories.
#
# scriptName - name of PBS file to create (including path)
# inputDirectory - full path to directory to process

def fastqc_PBS(scriptName, inputDirectory):

    # create file
    fo = open(scriptName, "w")
    
    # Write the file
    fo.write("#!/bin/bash\n")

    fo.write("#PBS -l walltime=06:00:00\n")   # wall time
    fo.write("#PBS -l ncpus=8,mem=3gb\n")     # PBS resources
    fo.write("#PBS -q uv2000\n")              # queue name

    # Command to find FASTQ files
    fo.write('FASTQ_FILES=`find ')
    fo.write(inputDirectory)
    fo.write(' -name "*.fastq.gz"` \n')

    # Call to run fastqc
    fo.write('fastqc --nogroup -t 8 -q $FASTQ_FILES \n')
    
    # close the PBS script
    fo.close()

