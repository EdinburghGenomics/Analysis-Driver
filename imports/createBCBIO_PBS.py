#!/opt/anaconda/bin/python

# Creates PBS scripts to run BCBIO
# INPUT:
#      1- i: Job number 
#      2- k: keys(sampleIds)
#      3- v: values -> tuple (Lane,Position,sampleName)

def bcbio_PBS(i, d, inputDirectory):

    from operator import itemgetter

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
    stringName ="#PBS -o bcbio"+`i`

    fo.write(stringName);
    fo.write("\n\n");

    # working directory
    fo.write("cd $PBS_O_WORKDIR\n\n");
    
    # base path to the BCL output
    # TODO: check the output contains "/Data/Intensities/BaseCalls"
    BASE_PATH= inputDirectory+"/Data/Intensities/BaseCalls"
    
    # get path to fastq pairs
    for k, v in sorted(d.items(), key=itemgetter(0)):
        #print 'id: %s, lanes: %s' % (k, v)
        FASTQ = BASE_PATH
    
        for j in xrange(0,len(v)):
            # check if sampleName != sampleId and if sampleName exists
            if v[0][2] != k and v[0][2] !='':
                FASTQ = BASE_PATH+"/"+str(v[0][2])+"/"+str(k)
                FASTQ1 =FASTQ+"/"+str(v[0][2])+"_S"+str(j+1)+"_L00"+str(v[j][0])+"_R"+str(j+1)+"_001"
#                print "[1]"+FASTQ1
                # sampleId = sampleName (must be != empty)
            elif v[0][2] == k and v[0][2] !='':
                FASTQ  = BASE_PATH+"/"+str(k)
                FASTQ1 = FASTQ+"/"+str(v[0][2])+"_S"+str(j+1)+"_L00"+str(v[j][0])+"_R"+str(j+1)+"_001"
 #               print "[2]"+FASTQ1
            else:
                FASTQ  = BASE_PATH+"/"+str(k)
                FASTQ1 = FASTQ + "/"+ str(k)+"_S"+str(j+1)+"_L00"+str(v[j][0])+"_R"+str(j+1)+"_001"
                
                
    fo.write("\n\n");
    
    # path to BCBIO
    BCBIO="/home/U008/lcebaman/bcbio/bin/bcbio_nextgen.py"
    PROJECT_PATH = "/home/U008/lcebaman/scripts/pbs"
    PROJECT_NAME = "project"+`i`
    PROJECT = PROJECT_PATH+"/"+PROJECT_NAME
    PROJECT_FLAGS="-w template gatk-variant"
    PROJECT_RUN = BCBIO +" "+ PROJECT_FLAGS +" "+ PROJECT_NAME +" "+ FASTQ1+".fastq.gz" +" "+ FASTQ2+".fastq.gz"
    
    # generate project to run bcbio
    fo.write(PROJECT_RUN);
    fo.write("\n\n");
    # bash command to run bcl2fastq
    BCBIO_RUN = BCBIO +" "+ PROJECT_NAME + "/config/"+PROJECT_NAME+".yaml" + " -n 16 " + \
        "--workdir" +" "+ PROJECT_NAME+"/work"
    fo.write(BCBIO_RUN);
    fo.write("\n");

    # close the PBS script
    fo.close();


def bcbio_loop(d, inputDirectory):

    # number of different PBS scripts
    n = len(d)
    
    BASE_PATH= inputDirectory+"/Data/Intensities/BaseCalls"

    for counter in xrange(0, n):
        # generate a PBS script per n
        bcbio_PBS(counter, d, inputDirectory)
        counter = counter + 1
    

#Unit Testing
# creates files from 0 to 2-1

