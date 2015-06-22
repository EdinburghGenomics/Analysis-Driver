def bcbio_pbs(sample_id, v, inputDirectory, project_name, sample_project):
    '''
    Writes out PBS bash scripts to run bcbio
    
    :param k: keys(sample ids)
    :param v: values -> tuple (Lane,Position,sample_name)
    3- inputDirectory: Path to input directory
    4- projectName: Name of the project
    5- sampleProject: Sample_Project
    '''

    # create a PBS script to run BCL2FASTQ
    lane = str(v[0])
    sample_name = str(v[1])
    pos = str(v[2])

    filename = project_name + '/pbs/runBCBIO_' + lane + '.pbs'
    fo = open(filename, 'w')

    fo.write('#!/bin/bash\n')

    fo.write('#PBS -l walltime=72:00:00\n')  # walltime needed
    fo.write('#PBS -l ncpus=8,mem=64gb\n')  # PBS resources
    fo.write('#PBS -N bcl2fastq\n')  # job name
    fo.write('#PBS -q uv2000\n')  # queue name
    fo.write('#PBS -j oe \n')  # input/output
    fo.write('#PBS -o ' + sample_name + '_' + lane + '\n\n')  # output file name

    fo.write('cd $PBS_O_WORKDIR\n\n\n')  # working directory

    # path to java
    fo.write('export JAVA_HOME=/home/U008/lcebaman/jdk1.7.0_76/\n')
    fo.write('export JAVA_BINDIR=/home/U008/lcebaman/jdk1.7.0_76/bin\n')

    # path to bcbio
    bcbio_home = '/home/U008/lcebaman/bcbio/bin'  # TODO: make this non-hard-coded
    bcbio_nextgen = bcbio_home + '/bcbio_nextgen.py'
    fastqc = bcbio_home + '/fastq'

    # base path to the BCL output
    # TODO: check the output contains '/Data/Intensities/BaseCalls'
    # base_path= inputDirectory +'/'+ sampleProject+'/'+sample_id  #+'/Data/Intensities/BaseCalls'

    base_path = '../Unaligned/' + sample_project + '/' + sample_id
    fastq1 = base_path + '/' + sample_name + '_S' + pos + '_L00' + lane + '_R1_001.fastq.gz'
    fastq2 = base_path + '/' + sample_name + '_S' + pos + '_L00' + lane + '_R2_001.fastq.gz'

    bcbio_gatk_template = bcbio_nextgen + ' -w template gatk-variant ' + sample_name + '_' + lane + ' ' + fastq1 + \
        ' ' + fastq2
    fo.write(bcbio_gatk_template + '\n\n')  # generate project to run bcbio

    # bash command to run bcl2fastq
    bcbio_run = bcbio_nextgen + ' ' + sample_name + '_' + lane + '/config/' + sample_name + '_' + lane + '.yaml' +\
        ' -n 16 --workdir ' + sample_name + '_' + lane + '/work\n\n'
    fo.write(bcbio_run)

    fo.write(fastqc + ' ' + '--nogroup -t 16' + ' -q ' + fastq1 + '-o ' + base_path + '\n')
    fo.write(fastqc + ' ' + '--nogroup -t 16' + ' -q ' + fastq2 + '-o ' + base_path + '\n\n')

    fo.close()


def bcbio_loop(d, input_dir, project_name, sample_project):
    """
    Creates a pbs script per pair of samples
    :param d:
    :param input_dir:
    :param project_name:
    :param sample_project:
    :return:
    """
    # number of different PBS scripts
    base_path = input_dir + '/Data/Intensities/BaseCalls'

    for k, v in d.items():
        for j in v:
            # generate a PBS script per n
            bcbio_pbs(k, j, input_dir, project_name, sample_project)
