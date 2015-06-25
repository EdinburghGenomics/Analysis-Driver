def fastqc_PBS(script_name, input_directory):
    """
    Write out a PBS Bash script to run fastqc on all fastq files in the specified directory
    :param script_name: Name of pbs file to create
    :param input_directory: full path to the dir to process
    :return:
    """

    fo = open(script_name, 'w')
    
    fo.write('#!/bin/bash\n')

    fo.write('#PBS -l walltime=06:00:00\n')  # wall time
    fo.write('#PBS -l ncpus=8,mem=3gb\n')  # PBS resources
    fo.write('#PBS -q uv2000\n')  # queue name

    # find FASTQ files
    fo.write('FASTQ_FILES=`find ' + input_directory + ' -name \'*.fastq.gz\'`\n')

    # run fastqc
    fo.write('fastqc --nogroup -t 8 -q $FASTQ_FILES\n')

    fo.close()
