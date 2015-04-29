#!/opt/anaconda/bin/python

# Creates PBS scripts to run BCBIO
# INPUT:
#      1- k: keys(sampleIds)
#      2- v: values -> tuple (Lane,Position,sampleName)
#      3- inputDirectory: Path to input directory
#      4- projectName: Name of the project

def bcbio_PBS( k, v, inputDirectory, projectName):

    # create a PBS script to run BCL2FASTQ
    sampleName = str(v[1])    
    lane = str(v[0])
    pos = str(v[2])

    filename=projectName+"/pbs/runBCBIO_"+lane+".pbs"
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
    stringName ="#PBS -o "+sampleName+"_"+lane

    fo.write(stringName);
    fo.write("\n\n");

    # working directory
    fo.write("cd $PBS_O_WORKDIR\n\n");
    
    # base path to the BCL output
    # TODO: check the output contains "/Data/Intensities/BaseCalls"
    BASE_PATH= inputDirectory+"/Data/Intensities/BaseCalls"
    
    # get path to fastq pairs
    FASTQ = BASE_PATH
    #    for j in xrange(0,len(v)):
        # check if sampleName != sampleId and if sampleName exists
    FASTQ = BASE_PATH+"/"+str(sampleName)+"/"+str(k)
    FASTQ1 =FASTQ+"/"+str(sampleName)+"_S"+str(pos)+"_L00"+str(lane)+"_R1_001"
    FASTQ2 =FASTQ+"/"+str(sampleName)+"_S"+str(pos)+"_L00"+str(lane)+"_R2_001"
                
    fo.write("\n\n");
    
    # path to BCBIO

    BCBIO="/home/U008/lcebaman/bcbio/bin/bcbio_nextgen.py"
    PROJECT_PATH = "/home/U008/lcebaman/scripts/pbs"
    PROJECT_NAME = "project"+lane
    PROJECT = PROJECT_PATH+"/"+PROJECT_NAME
    projFlags ="-w template gatk-variant"
    PROJECT_RUN = BCBIO +" "+ projFlags +" "+ sampleName +"_"+lane +" "+ FASTQ1+".fastq.gz" +" "+ FASTQ2+".fastq.gz"
     
    # generate project to run bcbio
    fo.write(PROJECT_RUN);
    fo.write("\n\n");
     # bash command to run bcl2fastq
    BCBIO_RUN = BCBIO +" "+ sampleName+"_"+lane + "/config/"+sampleName+"_"+lane+".yaml" + " -n 16 " + \
        "--workdir" +" "+ sampleName+"_"+lane +"/work"
    fo.write(BCBIO_RUN);
    fo.write("\n");
     
    # close the PBS script
    fo.close();


def bcbio_loop(d, inputDirectory, projectName):

    # number of different PBS scripts
    BASE_PATH= inputDirectory+"/Data/Intensities/BaseCalls"
 
    for k, v in d.items():
        for j in v:
            # generate a PBS script per n
            bcbio_PBS( k, j, inputDirectory, projectName)
              



