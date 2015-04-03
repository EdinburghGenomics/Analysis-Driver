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


def createBcl2fastq_PBS(mask):
    # create a PBS script to run BCL2FASTQ
    fo = open("pbs/run.pbs", "wb")

    fo.write("#!/bin/bash\n");
    # walltime needed
    fo.write( "#PBS -l walltime=00:10:00\n");

    # PBS resources
    fo.write("#PBS -l ncpus=4,mem=6gb\n");

    # jobname
    fo.write("#PBS -N bcl2fastq\n");

    # queue name
    fo.write("#PBS -q uv2000\n");

    # input/output
    fo.write("#PBS -j oe \n");

    # output file name
    fo.write("#PBS -o gccTest \n");

    # working directory
    fo.write("cd $PBS_O_WORKDIR\n");
    maskString = '--use-mask '
    maskString += generateMask(mask)

    # bash command to run bcl2fastq
    fo.write("bscl2fastq -r blabla -o bloblo ");
    fo.write(maskString);
    

    # close the PBS script
    fo.close()


#Unit Testing

# generate fake masklist
masklist=['128','Y','8','N','128','Y']
# call to generate string
maskString = generateMask(masklist)
print maskString


createBcl2fastq_PBS(masklist)    
