#!/opt/anaconda/bin/python

# INPUT:
#      1- i: Project and job number (total of 36-1 job numbers)
#      2- fastq1: First file of pair fastq files
#      3- fastq2: Second pair of pair fastq files
def bcbio_PBS(i, k, v):
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
    

    BASE_PATH= inputDirectory+"/Data/Intensities/BaseCalls"
    FASTQ_list=[]
    for k, v in sorted(d.items(), key=itemgetter(0)):

        # loop over the element values of the dictionary
        for i in xrange(0,len(v)):                   
        if v[0][2] != '':
            FASTQ=  BASE_PATH+"/"+str(v[0][2])+"/"+str(k)+"/"
            FASTQ =FASTQ+"/"+str(v[0][2])+"_S"+str(i+1)+"_L00"+str(v[i][0])+"_R"+str(i+1)+"_001"
            FASTQ2 =FASTQ+"/"+str(v[0][2])+"_S"+str(i+1)+"_L00"+str(v[i][0])+"_R"+str(i+1)+"_001"
            FASTQ_2 = FASTQ+str(v[0][2])+"_S"+str(v[1][1])+"_L00"+str(v[1][0])+"_R2"+"_001"
        else:
            FASTQ_1 = BASE_PATH+"/"+str(k)+"/"+str(k)+"_S"+str(v[0][1])+"_L00"+str(v[0][0])+"_R1"+"_001"
            FASTQ_2 = BASE_PATH+"/"+str(k)+"/"+str(k)+"_S"+str(v[1][1])+"_L00"+str(v[1][0])+"_R2"+"_001"


    # untar fastq files
    tarcom ="gunzip"+""+fastq1+"/"+"*.gz"



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


def bcbio_loop(d, inputDirectory):

    # number of different PBS scripts
    n = len(d)
    
    BASE_PATH= inputDirectory+"/Data/Intensities/BaseCalls"

    for counter in xrange(0, n):
        # generate a PBS script per n
        bcbio_PBS(counter, k , v )
        counter = counter + 1
    

#Unit Testing
# creates files from 0 to 2-1

