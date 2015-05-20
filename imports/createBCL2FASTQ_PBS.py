#!/opt/anaconda/bin/python

#Generate the mask as it is understood by BCL2FASTQ
# INPUT: A mask list with made of triplets
# OUTPUT: Mask string ready for BCL2FASTQ
#
# The string is built in reverse, so we can easily remove the last base

def generateMask(mask):
    masklen = len(mask)

    chain=[]

    # mask is a list of triplets, so we reverse step in threes
    currentReadNumber=0
    for i in xrange(masklen-3,-1,-3):

        readNumber=mask[i]
        numCycles=mask[i+1]
        isIndexedRead=mask[i+2]

        # build the mask string in reverse
        endChar=''
        if (readNumber != currentReadNumber):
            # New read
            currentReadNumber=readNumber
            numCycles=str(int(numCycles)-1)
            endChar='n'
        
        maskPart=numCycles+endChar
        
        if(isIndexedRead.lower() == 'n'):
            maskPart='y'+maskPart
        else:
            maskPart='i'+maskPart

        chain.append(maskPart)

    return ','.join(chain[::-1])

#Create the PBS script responsible for running BCL2FASTQ
#INPUT: A mask list

def bcl2fastq_PBS(mask, scriptName, projectName, inputDirectory):
    # create a PBS script to run BCL2FASTQ
    PBSname = projectName+"/pbs/"+scriptName
    fo = open(PBSname, "wb")

    fo.write("#!/bin/bash\n");
    # walltime needed
    fo.write( "#PBS -l walltime=24:00:00\n");

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
#masklist=['1','128','Y','2','8','N','3','128','Y']
# call to generate string
#maskString = generateMask(masklist)
#print maskString
# test PBS script generation for BCL2FASTQ
#createBcl2fastq_PBS(masklist)    
