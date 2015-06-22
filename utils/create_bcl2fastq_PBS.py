def generate_mask(mask):
    """
    Generate the mask as it is understood by BCL2FASTQ
    INPUT: A mask list with made of triplets
    OUTPUT: Mask string ready for BCL2FASTQ

    The string is built in reverse, so we can easily remove the last base
    :param mask: a list of triplets:
    :return:
    """

    chain = []

    # mask is a list of triplets, so we reverse step in threes
    current_read_number = 0
    for i in range(len(mask) - 3, -1, -3):

        read_number = mask[i]
        num_cycles = mask[i + 1]
        is_indexed_read = mask[i + 2]

        # build the mask string in reverse
        end_char = ''
        if read_number != current_read_number:
            # new read
            current_read_number = read_number
            num_cycles = str(int(num_cycles) - 1)
            end_char = 'n'
        
        mask_part = num_cycles+end_char
        
        if is_indexed_read.lower() == 'n':
            mask_part = 'y' + mask_part
        else:
            mask_part = 'i' + mask_part

        chain.append(mask_part)

    return ','.join(chain[::-1])


def bcl2fastq_PBS(mask, script_name, project_name, input_dir):
    """
    Create the PBS script responsible for running BCL2FASTQ
    :param mask: A mask list
    :param script_name:
    :param project_name:
    :param input_dir:
    :return:
    """

    PBS_name = project_name + '/pbs/' + script_name
    fo = open(PBS_name, 'w')  # create a PBS script to run BCL2FASTQ

    fo.write('#!/bin/bash\n')
    # walltime needed
    fo.write('#PBS -l walltime=24:00:00\n')

    # PBS resources
    fo.write('#PBS -l ncpus=12,mem=24gb\n')
    # jobname
    fo.write('#PBS -N bcl2fastq\n')

    # queue name
    fo.write('#PBS -q uv2000\n')

    # input/output
    fo.write('#PBS -j oe\n')

    # output file name
    fo.write('#PBS -o bcl2fastq.out\n\n')

    # working directory
    fo.write('cd $PBS_O_WORKDIR\n\n')

    # TODO: understand the mask issue
    mask_string = '--use-mask ' + generate_mask(mask)
    input_opt = '-R ' + input_dir
    
    # TODO: include the right bcl2fastq command
    # bash command to run bcl2fastq 
    # fo.write('dd if=/dev/zero of=/scratch/U008/lcebaman/test.bla  bs=32768 count=100000')
    fo.write('/scratch/U008/edingen/bin/bcl2fastq_2_16/bin/bcl2fastq ' + input_opt + ' -o ../Unaligned ')

    fo.write('--sample-sheet ' + input_dir + '/SampleSheet.csv ')
    # fo.write(mask_string)

    fo.close()


# Unit Testing
if __name__ == '__main__':
    # generate fake masklist
    masklist = ['1', '128', 'Y', '2', '8', 'N', '3', '128', 'Y']
    # call to generate string
    mask_string = generate_mask(masklist)
    print(mask_string)
    # test PBS script generation for BCL2FASTQ
    # createBcl2fastq_PBS(masklist)
