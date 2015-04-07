#!/opt/anaconda/bin/python

# INPUT:
#      1- i: Project and job number (total of 36-1 job numbers)
#      2- fastq1: First file of pair fastq files
#      3- fastq2: Second pair of pair fastq files
def bcbio_PBS(i, fastq1, fastq2):
    # create a PBS script to run BCL2FASTQ
    filename="pbs/runBCBIO"+`i`+".pbs"
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
    stringName ="#PBS -o bcbio"+`i`

    fo.write(stringName);
    fo.write("\n\n");

    # working directory
    fo.write("cd $PBS_O_WORKDIR\n\n");

    # path to BCBIO
    BCBIO="/home/U008/lcebaman/bcbio/bin/bcbio_nextgen.py"
    PROJECT_PATH = "/home/U008/lcebaman/scripts/pbs"
    PROJECT_NAME = "project"+`i`
    PROJECT = PROJECT_PATH+"/"+PROJECT_NAME
    PROJECT_FLAGS="-w template gatk-variant"
    PROJECT_RUN = BCBIO +" "+ PROJECT_FLAGS +" "+ PROJECT_NAME +" "+ fastq1 +" "+ fastq2
    
    # generate project to run bcbio
    fo.write(PROJECT_RUN);
    fo.write("\n\n");
    # bash command to run bcl2fastq
    BCBIO_RUN = BCBIO +" "+ PROJECT_NAME + "/config/"+PROJECT_NAME+".yaml" + " -n 16 " + "--workdir" +" "+ PROJECT_NAME+"/work"
    fo.write(BCBIO_RUN);
    fo.write("\n");

    # close the PBS script
    fo.close()

# TODO: We will probably need some more input parameters here since we have to use the fastq files to create the
#       the BCBIO projects. It is not clear the number of output fastq files that we will have and it supposes to
#       be the same number as the number of runs (36), otherwise we will be submmitting the same job more than once
def bcbio_loop(n):
    
    for i in xrange(0,n):
        # TODO: pass info about the fastq files
        FASTQ_PATH="/scratch/U008/lcebaman/fastq/1"
        FASTQ1 = FASTQ_PATH+"/1--NA12878-OD1_1.fastq"
        FASTQ2 = FASTQ_PATH+"/1--NA12878-OD1_2.fastq"
        #generate a PBS script per n
        bcbio_PBS(i, FASTQ1, FASTQ2 )
    

#Unit Testing
# creates files from 0 to 2-1

