def bcbio_PBS(k, v, input_directory, project_name, sample_project):
    """
    Creates PBS scripts to run bcbio
    :param k: keys (sampleIDs)
    :param v: values -> tuple(Lane, Position, sample_name)
    :param input_directory
    :param project_name
    :param sample_project 
    :return:
    """

    # create a PBS script to run BCL2FASTQ
    sample_name = str(v[1])    
    lane = str(v[0])
    pos = str(v[2])
    sample_id = str(k)

    file_name = project_name + '/pbs/runBCBIO_'+ lane + '.pbs'
    fo = open(file_name, 'w')

    fo.write('#!/bin/bash\n')
    fo.write('#PBS -l walltime=72:00:00\n')  # walltime needed
    fo.write('#PBS -l ncpus=8,mem=64gb\n')  # PBS resources
    fo.write('#PBS -N bcl2fastq\n')  # job name
    fo.write('#PBS -q uv2000\n')  # queue name
    fo.write('#PBS -j oe \n')  # input/output
    fo.write('#PBS -o ' + sample_name + '_' + lane)  # output file name
    fo.write('\n\n')

    fo.write('cd $PBS_O_WORKDIR\n\n')  # working directory
    
    # base path to the BCL output
    # TODO: check the output contains '/Data/Intensities/BaseCalls'
    # base_path = inputDirectory + '/' + sampleProject + '/' + sample_id  # + '/Data/Intensities/BaseCalls'
    base_path = '../Unaligned' + '/' + sample_project + '/' + sample_id
    fastq1 = base_path + '/' + sample_name + '_S' + pos + '_L00' + lane + '_R1_001.fastq.gz'
    fastq2 = base_path + '/' + sample_name + '_S' + pos + '_L00' + lane + '_R2_001.fastq.gz'

    fo.write('export JAVA_HOME=/home/U008/lcebaman/jdk1.7.0_76/\n')  # paths to java
    fo.write('export JAVA_BINDIR=/home/U008/lcebaman/jdk1.7.0_76/bin\n')

    # path to bcbio
    bcbio_home = '/home/U008/lcebaman/bcbio/bin'
    bcbio = bcbio_home + '/bcbio_nextgen.py'
    fastqc = bcbio_home + '/fastq'

    proj_flags = '-w template gatk-variant'
    bcbio_variant_call = bcbio + ' ' + proj_flags + ' ' + sample_name + '_' + lane + ' ' + fastq1 + ' ' + fastq2

    fo.write(bcbio_variant_call + '\n\n')  # generate project to run bcbio
    # bash command to run bcl2fastq
    bcbio_run = bcbio + ' ' + sample_name + '_' + lane + '/config/' + sample_name + '_' + lane + '.yaml' + ' -n 16 ' + \
        '--workdir' + ' ' + sample_name + '_' + lane + '/work'

    fo.write(bcbio_run + '\n\n')

    fo.write(fastqc + ' ' + '--nogroup -t 16' + ' -q ' + fastq1 + '-o ' + base_path + '\n')
    fo.write(fastqc + ' ' + '--nogroup -t 16' + ' -q ' + fastq2 + '-o ' + base_path + '\n')

    fo.close()


def bcbio_loop(d, input_directory, project_name, sample_project):
    """
    Create a PBS script per pair of samples
    :param d:
    :param input_directory:
    :param project_name:
    :param sample_project:
    :return:
    """
 
    for k, v in d.items():
        for j in v:
            # generate a PBS script per n
            bcbio_PBS(k, j, input_directory, project_name, sample_project)
