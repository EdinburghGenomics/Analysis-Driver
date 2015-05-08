#!/opt/anaconda/bin/python

# Creates PBS scripts to run BCBIO
# INPUT:
#      1- k: keys(sampleIds)
#      2- v: values -> tuple (Lane,Position,sampleName)
#      3- inputDirectory: Path to input directory
#      4- projectName: Name of the project
#      5- sampleProject: Sample_Project 
def bcbio_PBS( k, v, inputDirectory, projectName, sampleProject):

    # create a PBS script to run BCL2FASTQ
    sampleName = str(v[1])    
    lane = str(v[0])
    pos = str(v[2])
    sampleID=str(k)

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
    #BASE_PATH= inputDirectory +"/"+ sampleProject+"/"+sampleID  #+"/Data/Intensities/BaseCalls"
    BASE_PATH="../Unaligned"+"/"+ sampleProject+"/"+sampleID
    FASTQ1 = BASE_PATH +"/"+str(sampleName)+"_S"+str(pos)+"_L00"+str(lane)+"_R1_001.fastq.gz"
    FASTQ2 = BASE_PATH +"/"+str(sampleName)+"_S"+str(pos)+"_L00"+str(lane)+"_R2_001.fastq.gz"
                
    fo.write("\n\n");

    #path to java
    fo.write("export JAVA_HOME=/home/U008/lcebaman/jdk1.7.0_76/");
    fo.write("export JAVA_BINDIR=/home/U008/lcebaman/jdk1.7.0_76/bin");

    # path to BCBIO
    BCBIO_HOME="/home/U008/lcebaman/bcbio/bin"
    BCBIO=BCBIO_HOME+"/bcbio_nextgen.py"
    FASTQC=BCBIO_HOME+"/fastq"

    projFlags ="-w template gatk-variant"
    PROJECT_RUN = BCBIO +" "+ projFlags +" "+ sampleName +"_"+lane +" "+ FASTQ1 +" "+ FASTQ2
     
    # generate project to run bcbio
    fo.write(PROJECT_RUN);
    fo.write("\n\n");
     # bash command to run bcl2fastq
    BCBIO_RUN = BCBIO +" "+ sampleName+"_"+lane + "/config/"+sampleName+"_"+lane+".yaml" + " -n 16 " + \
        "--workdir" +" "+ sampleName+"_"+lane +"/work"

    fo.write(BCBIO_RUN);
    
    fo.write("\n\n");
    fastqArgs = FASTQC + " " + "--nogroup -t 16"+ " -q "+ FASTQ1 + "-o " + BASE_PATH
    fo.write(fastqArgs);
    fo.write("\n\n");
    fastqArgs = FASTQC + " " + "--nogroup -t 16"+ " -q "+ FASTQ2 + "-o " + BASE_PATH
    fo.write(fastqArgs);

    fo.write("\n\n");
     
    # close the PBS script
    fo.close();

# Creates a pbs script per pair of samples
def bcbio_loop(d, inputDirectory, projectName, sampleProject):

    # number of different PBS scripts
    BASE_PATH= inputDirectory+"/Data/Intensities/BaseCalls"
 
    for k, v in d.items():
        for j in v:
            # generate a PBS script per n
            bcbio_PBS( k, j, inputDirectory, projectName, sampleProject)
              



