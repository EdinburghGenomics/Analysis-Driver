def bcbio_PBS(k, v, input_dir, pbs_name, sample_project):
    """
    Writes out a PBS Bash script to run bcbio
    :param k: key - sample ID
    :param v: values -> tuple(lane, position, sample_name)
    :param input_dir
    :param pbs_name
    :param sample_project
    :return:
    """

    # create a PBS script to run BCL2FASTQ
    sample_name = str(v[1])
    if v[0]:
        lane = '_L00' + v[0]
    else:
        lane = ''

    if v[2]:
        pos = '_S' + v[0]
    else:
        pos = ''

    sample_id = k

    f = open(pbs_name, 'w')

    f.write('#!/bin/bash\n')

    f.write('#PBS -l walltime=72:00:00\n')  # walltime needed
    f.write('#PBS -l ncpus=8,mem=64gb\n')  # PBS resources
    f.write('#PBS -N bcbio\n')  # jobname
    f.write('#PBS -q uv2000\n')  # queue name
    f.write('#PBS -j oe\n')  # input/output
    f.write('#PBS -o ' + sample_name + lane)  # output file name
    f.write('\n\n')

    f.write('cd $PBS_O_WORKDIR\n\n')  # working directory
                
    f.write('\n\n')

    # paths to java
    f.write('export JAVA_HOME=/home/U008/lcebaman/jdk1.7.0_76/\n')
    f.write('export JAVA_BINDIR=/home/U008/lcebaman/jdk1.7.0_76/bin\n')

    # path to bcbio
    bcbio_home = '/home/U008/lcebaman/bcbio/bin'
    bcbio = bcbio_home + '/bcbio_nextgen.py'
    fastqc = bcbio_home + '/fastq'

    # base path to the BCL output
    # TODO: check the output contains '/Data/Intensities/BaseCalls'
    # base_path = inputDirectory + '/' + sampleProject + '/' + sample_id  # +'/Data/Intensities/BaseCalls'
    base_path = '../fastq/' + sample_project + '/' + sample_id
    fastq1 = base_path + '/' + sample_name + pos + lane + '_R1_001.fastq.gz'
    fastq2 = base_path + '/' + sample_name + pos + lane + '_R2_001.fastq.gz'

    bcbio_template = bcbio + ' -w template gatk-variant ' + sample_name + '_' + lane + ' ' + fastq1 + ' ' + fastq2

    # bash command to set up bcbio project
    f.write(bcbio_template)
    f.write('\n\n')

    # bash command to run bcbio
    bcbio_run = bcbio + ' ' + sample_name + '_' + lane + '/config/' + sample_name + '_' + lane + '.yaml' + ' -n 16 ' + \
        '--workdir ' + sample_name + '_' + lane + '/work'

    f.write(bcbio_run)

    f.write('\n\n')
    f.write(fastqc + ' --nogroup -t 16' + ' -q ' + fastq1 + ' -o ' + base_path + '\n')
    f.write(fastqc + ' ' + '--nogroup -t 16'+ ' -q '+ fastq2 + '-o ' + base_path + '\n')

    f.close()


# Creates a pbs script per pair of samples
def bcbio_loop(d, input_dir, pbs_name, sample_project):
    for k, v in d.items():
        for j in v:
            # generate a PBS script per n
            bcbio_PBS(k, j, input_dir, pbs_name, sample_project)


