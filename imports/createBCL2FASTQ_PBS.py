#!/opt/anaconda/bin/python

#Generate the mask as it is understood by BCL2FASTQ
# INPUT: A mask list with made of pairs
# OUTPUT: Mask string ready for BCL2FASTQ

def generateMask(mask):
    masklen = len(mask)

    chain=''
    # mask is a list of pairs, so we just need length/2
    for i in xrange(0,(masklen/2)):
        j =i*2
        # need to add a comma if not the last
        if (j != masklen-2):
            chain += mask[j]+mask[j+1]+','
        else:
            chain += mask[j]+mask[j+1]

    return chain

#Create the PBS script responsible for running BCL2FASTQ
#INPUT: A mask list

def bcl2fastq_PBS(mask, scriptName, projectName, inputDirectory):
    # create a PBS script to run BCL2FASTQ
    PBSname = projectName+"/pbs/"+scriptName
    fo = open(PBSname, "wb")

    fo.write("#!/bin/bash\n");
    # walltime needed
    fo.write( "#PBS -l walltime=03:00:00\n");

    # PBS resources
    fo.write("#PBS -l ncpus=12,mem=24gb\n");

    # jobname
    fo.write("#PBS -N bcl2fastq\n");

    # queue name
    fo.write("#PBS -q uv2000\n");

    # input/output
    fo.write("#PBS -j oe \n");

    # output file name
    fo.write("#PBS -o bcl2fastq.out \n\n");

    # working directory
    fo.write("cd $PBS_O_WORKDIR\n\n");

    #TODO: understand the mask issue
    maskString = '--use-mask '
    maskString += generateMask(mask) 
    inputOpt = "-R " + inputDirectory
    
    # TODO: include the right bcl2fast command
    # bash command to run bcl2fastq 
    #fo.write("dd if=/dev/zero of=/scratch/U008/lcebaman/test.bla  bs=32768 count=100000");
    BCLPATH = '/scratch/U008/edingen/bin/bcl2fastq_2_16/bin/bcl2fastq'
    fo.write(BCLPATH);
    fo.write(" ");
    fo.write(inputOpt);
    fo.write(" ");
    outputOpt = "-o ../Unaligned"
    fo.write(outputOpt);
    fo.write(" ");
    sampleSheetPath=inputDirectory+"/SampleSheet.csv"
    sampleOpt ="--sample-sheet "+ sampleSheetPath
    fo.write(sampleOpt);
    fo.write(" ");
    #fo.write(maskString);

    
    # close the PBS script
    fo.close()


#Unit Testing

# generate fake masklist
#masklist=['128','Y','8','N','128','Y']
# call to generate string
#maskString = generateMask(masklist)
#print maskString
# test PBS script generation for BCL2FASTQ
#createBcl2fastq_PBS(masklist)    
