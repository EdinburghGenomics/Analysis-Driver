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

def bcl2fastq_PBS(mask, pbsName, inputDirectory, fastqPath):
    # create a PBS script to run BCL2FASTQ
    fo = open(pbsName, "w")

    # script header
    fo.write("#!/bin/bash\n");

    # walltime needed
    fo.write( "#PBS -l walltime=24:00:00\n");

    # PBS resources
    fo.write("#PBS -l ncpus=12,mem=24gb\n");

    # queue name
    fo.write("#PBS -q uv2000\n");

    # construct command and arguments
    fo.write('bcl2fastq ')
    fo.write(' --runfolder-dir '  + inputDirectory)
    fo.write(' --output-dir '     + fastqPath)
    fo.write(' --sample-sheet '   + inputDirectory+'SampleSheet.csv')
    fo.write(' --use-bases-mask ' + generateMask(mask))

    # close the PBS script
    fo.write("\n")
    fo.close()


#Unit Testing

# generate fake masklist
#masklist=['1','128','Y','2','8','N','3','128','Y']
# call to generate string
#maskString = generateMask(masklist)
#print maskString
# test PBS script generation for BCL2FASTQ
#createBcl2fastq_PBS(masklist)    
