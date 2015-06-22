#!/opt/anaconda/bin/python

# This program expects an input directory. From there, we will run BCL2FASTQ and BCBIO
import sys
import os
import logging
from datetime import datetime

from utils import run_info_parser
from utils import sample_sheet_parser
from utils import create_bcl2fastq_PBS
from utils import create_bcbio_PBS
from utils import qsub_dependents
from utils import args
from utils import make_project


if __name__ == '__main__':

    # Open logging file and configure it with date and time stamps
    log_file_name = 'log_' + datetime.now().strftime('%Y%m%d-%H%M%S') + '.log'
    logging.basicConfig(filename=log_file_name, format='%(asctime)s %(message)s', datefmt='[%d/%m/%Y-%H:%M:%S]')

    # parse the input directory
    input_path = args.get_file_name(sys.argv[1:])
    
    logging.info('Reading bcl data from ' + input_path)

    # create project directory
    project_name = make_project.get_dirname(input_path)
    work_dir = project_name + '_work'
    make_project.make_project(work_dir)

    # Read RunInfo.xml
    logging.info('Reading the mask from ' + input_path)
    mask = run_info_parser.get_mask(input_path)
    
    # Read SampleSheet.csv
    logging.info('Reading SampleSheet from %s ' + input_path)
    num_lanes, sheet_dict = sample_sheet_parser.read_sample_sheet(input_path)
    logging.info('Lanes: ' + num_lanes)
    # get sample_project
    sample_project = sample_sheet_parser.get_sample_project(sheet_dict)

    # Create BCL2FASTQ PBS script
    logging.info('Creating BCL2FASTQ pbs script')
    psb_name = 'BCL_' + project_name + '.pbs'
    logging.info('name: ' + psb_name)

    create_bcl2fastq_PBS.bcl2fastq_PBS(mask, psb_name, work_dir, input_path)

    # create BCBIO PBS scripts
    logging.info('Creating BCBIO PBS scripts')
    create_bcbio_PBS.bcbio_loop(sheet_dict, input_path, work_dir, sample_project)
    
    os.chdir(work_dir + '/pbs')
    
    # submit bcl2fastq 
    # create a list with the name of the PBS script
    args = [psb_name]
    
    # submit the BCL2FASTQ script to batch scheduler
    logging.info('Submitting BCL2FASTQ_PBS')
    bcl2fastq_job_id = qsub_dependents.qsub_dependents(args)
    
    # submit the BCBIO scripts once BCL2FASTQ has finished
    logging.info('Submitting BCBIO_PBS')

    job_ids = []
    # submit set of BCBIO jobs. A job per sampleId included in the SampleSheet
    for i in range(1, num_lanes):
        script_name = 'runBCBIO_' + str(i) + '.pbs'
        args = [script_name]
        # store the jobIds in a list. They will not get executed until BCL2FASTQ has finished
        job_ids.append(qsub_dependents.qsub_dependents(args, bcl2fastq_job_id))
         
    logging.info('Submitted BCBIO_PBS')
