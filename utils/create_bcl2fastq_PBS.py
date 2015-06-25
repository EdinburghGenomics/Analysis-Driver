def generate_mask(mask):
    """
    Generate a mask understood by BCL2FASTQ. The string is built in reverse to easily remove the last base
    :param mask: A list of triplets
    :return: A mask string ready for bcl2fastq
    """
    chain = []

    # mask is a list of triplets, so we reverse step in threes
    current_read_number=0
    for i in range(len(mask) - 3, -1, -3):
        read_number = mask[i]
        num_cycles = mask[i + 1]
        is_indexed_read = mask[i + 2]

        # build the mask string in reverse
        end_char = ''
        if read_number != current_read_number:
            # new read
            current_read_number = read_number
            num_cycles = str(int(num_cycles)-1)
            end_char = 'n'
        
        mask_part = num_cycles + end_char
        
        if is_indexed_read.lower() == 'n':
            mask_part = 'y' + mask_part
        else:
            mask_part = 'i' + mask_part

        chain.append(mask_part)

    return ','.join(chain[::-1])


def bcl2fastq_PBS(mask, script_name, project_name, input_dir):
    """
    Create a PBS script o run bcl2fastq
    :param mask: A mask list
    :param script_name:
    :param project_name:
    :param input_dir:
    :return:
    """
    # create a PBS script to run BCL2FASTQ
    pbs_name = project_name + '/pbs/' + script_name
    f = open(pbs_name, 'w')

    f.write('#!/bin/bash\n')

    f.write('#PBS -l walltime=24:00:00\n')  # walltime needed
    f.write('#PBS -l ncpus=12,mem=24gb\n')  # PBS resources
    f.write('#PBS -N bcl2fastq\n')  # job name
    f.write('#PBS -q uv2000\n')  # queue name
    f.write('#PBS -j oe\n')  # input/output
    f.write('#PBS -o bcl2fastq.out\n\n')  # output file name

    f.write('cd $PBS_O_WORKDIR\n\n')  # working directory

    input_opt = '-R ' + input_dir
    
    # TODO: include the right bcl2fast command
    # bash command to run bcl2fastq 
    # f.write('dd if=/dev/zero of=/scratch/U008/lcebaman/test.bla  bs=32768 count=100000\n')
    bcl_path = '/scratch/U008/edingen/bin/bcl2fastq_2_16/bin/bcl2fastq'
    f.write(bcl_path + ' ' + input_opt + ' -o ../Unaligned --sample-sheet ' + input_dir + '/SampleSheet.csv ')
    # TODO: understand the mask issue
    # mask_string = '--use-mask '
    # mask_string += generate_mask(mask)
    # f.write(mask_string)

    f.close()


# Unit test
if __name__ == '__main__':
    # generate fake masklist
    masklist=['1','128','Y','2','8','N','3','128','Y']
    # call to generate string
    maskString = generate_mask(masklist)
    print(maskString)
    # test PBS script generation for BCL2FASTQ
    bcl2fastq_PBS(masklist)
