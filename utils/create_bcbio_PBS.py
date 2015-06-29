import os

def bcbio_PBS(pbs_name, bcbio_run_folder, run_id, lane, sample_name='Unassigned_S0'):
    """
    Writes out a PBS Bash script to run bcbio
    :param k: key - sample ID
    :param v: values -> tuple(lane, sample_name, position)
    :param input_dir
    :param pbs_name
    :param sample_project
    :return:
    """

    bcbio_run_folder += run_id + '_L00' + lane + '/'

    if not os.path.exists('../jobs/' + run_id + '/bcbio'):
        os.mkdir('../jobs/' + run_id + '/bcbio')
    pbs_name += '_L00' + lane + '.pbs'
    print('Opening ' + pbs_name)
    # create a PBS script to run BCL2FASTQ
    f = open(pbs_name, 'w')

    f.write('#!/bin/bash\n')

    f.write('#PBS -l walltime=72:00:00\n')  # walltime needed
    f.write('#PBS -l ncpus=8,mem=64gb\n')  # PBS resources
#    f.write('#PBS -N bcbio\n')  # jobname
    f.write('#PBS -q uv2000\n')  # queue name
    f.write('#PBS -j oe\n')  # input/output
    f.write('#PBS -o ' + bcbio_run_folder + 'out.txt')  # output file name
    f.write('\n\n')

    f.write('cd $PBS_O_WORKDIR\n\n')  # working directory
                
    # paths to java
    f.write('export JAVA_HOME=/home/U008/lcebaman/jdk1.7.0_76/\n')
    f.write('export JAVA_BINDIR=/home/U008/lcebaman/jdk1.7.0_76/bin\n\n')

    # path to bcbio
    bcbio_home = '/home/U008/lcebaman/bcbio/bin'
    bcbio = bcbio_home + '/bcbio_nextgen.py'
    fastqc = bcbio_home + '/fastqc'

    # Path to the input fastqs
    # base_path = inputDirectory + '/' + sampleProject + '/' + sample_id  # +'/Data/Intensities/BaseCalls'
    base_path = '../fastq/' + run_id
    fastq1 = base_path + '/' + sample_name + '_L00' + lane + '_R1_001.fastq.gz'
    fastq2 = base_path + '/' + sample_name + '_L00' + lane + '_R2_001.fastq.gz'

    bcbio_template = ' '.join(
        [
            bcbio,
            '-w', 'template', 'gatk-variant',
            bcbio_run_folder,
            fastq1,
            fastq2
        ]
    )

    # bash command to set up bcbio project
    f.write(bcbio_template)
    f.write('\n\n')

    print('Sample name: ' + sample_name)
    # bash command to run bcbio
    bcbio_run = ' '.join(
        [
            bcbio,
            bcbio_run_folder + 'config/' + sample_name + '_' + lane + '.yaml',
            '-n', '16',
            '--workdir', bcbio_run_folder + 'work'
        ]
    )

    f.write(bcbio_run)

    f.write('\n\n')
    # f.write(fastqc + ' --nogroup -t 16' + ' -q ' + fastq1 + ' -o ' + base_path + '\n')
    # f.write(fastqc + ' ' + '--nogroup -t 16'+ ' -q '+ fastq2 + '-o ' + base_path + '\n')

    f.close()


