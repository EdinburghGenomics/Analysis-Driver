#!/opt/anaconda/bin/python

# 
def createBCBIO_PBS(i):
    # create a PBS script to run BCL2FASTQ
    filename="pbs/runBCBIO_"+`i`+".pbs"
    fo = open(filename, "wb")

    fo.write("#!/bin/bash\n");
    # walltime needed
    fo.write( "#PBS -l walltime=72:00:00\n");

    # PBS resources
    fo.write("#PBS -l ncpus=8,mem=64gb\n");

    # jobname
    fo.write("#PBS -N bcl2fastq\n");

    # queue name
    fo.write("#PBS -q uv2000\n");

    # input/output
    fo.write("#PBS -j oe \n");

    # output file name
    stringName ="#PBS -o bcbio_"+`i`

    fo.write(stringName);
    fo.write("\n");

    # working directory
    fo.write("cd $PBS_O_WORKDIR\n");

    # bash command to run bcl2fastq
    fo.write("bcbio_nextq ......");

    # close the PBS script
    fo.close()


#Unit Testing
# creates files from 0 to 2-1
for i in xrange(0,2):
    createBCBIO_PBS(i)

