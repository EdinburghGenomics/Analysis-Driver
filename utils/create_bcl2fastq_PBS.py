def generate_mask(mask):
    """
    Generate a mask readable by bcl2fastq. The string is built in reverse to easily remove the last base
    :param mask: A list of triplets
    :return: A string readable by bcl2fastq
    """

    chain = []

    # mask is a list of triplets, so we reverse step in threes
    current_read_number=0
    for i in range(len(mask)-3, -1, -3):

        read_number = mask[i]
        num_cycles = mask[i+1]
        is_indexed_read = mask[i+2]

        # build the mask string in reverse
        end_char = ''
        if read_number != current_read_number:
            # New read
            current_read_number = read_number
            num_cycles = str(int(num_cycles)-1)
            end_char = 'n'
        
        mask_part = num_cycles+end_char
        
        if is_indexed_read.lower() == 'n':
            mask_part = 'y' + mask_part
        else:
            mask_part = 'i' + mask_part

        chain.append(mask_part)

    return ','.join(chain[::-1])


def bcl2fastq_PBS(mask, pbs_name, inputDirectory, fastq_path):
    """
    Write out a PBS Bash script for running bcl2fastq
    :param mask: Mask list
    :param pbs_name:
    :param inputDirectory:
    :param fastq_path:
    :return:
    """
    # create a PBS script to run BCL2FASTQ
    f = open(pbs_name, 'w')

    f.write('#!/bin/bash\n')  # script header
    f.write('#PBS -l walltime=24:00:00\n')  # walltime needed
    f.write('#PBS -l ncpus=12,mem=24gb\n')  # PBS resources
    f.write('#PBS -q uv2000\n')  # queue name

    # construct command and arguments
    f.write(
        'bcl2fastq --runfolder-dir ' + inputDirectory + ' --output-dir ' + fastq_path + ' --sample-sheet ' +
        inputDirectory + 'SampleSheet.csv --use-bases-mask ' + generate_mask(mask) + '\n'
    )

    f.close()


# unit test
if __name__ == '__main__':
    # generate fake masklist
    masklist=['1','128','Y','2','8','N','3','128','Y']
    # call to generate string
    maskString = generate_mask(masklist)
    print(maskString)
    # test PBS script generation for BCL2FASTQ
    bcl2fastq_PBS(masklist)
