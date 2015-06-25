#!/opt/anaconda/bin/python

# This program expects an input directory. From there, we will run BCL2FASTQ and BCBIO
import sys
import os
import logging
from datetime import datetime

from utils import run_info_parsing
from utils import sample_sheet_parsing
from utils import create_bcl2fastq_PBS
from utils import create_bcbio_PBS
from utils import qsub_dependents
from utils import args
from utils import make_project


if __name__ == '__main__':

    log_file = 'log_' + datetime.now().strftime('%Y%m%d-%H%M%S') + '.log'
    logging.basicConfig(
        filename=log_file,
        level=logging.DEBUG,
        format='%(asctime)s %(message)s',
        datefmt='[%d/%m/%Y-%H:%M:%S]'
    )
    logger = logging.getLogger('AnalysisDriver')

    # parse the input directory
    input_path = args.parse_args(sys.argv[1:])
    
    logger.info('Reading bcl data from %s ', input_path)

    # create project directory
    project_name = make_project.get_dir_name(input_path)
    work_dir = project_name + '_work'
    make_project.make_project(work_dir)

    logger.info('Reading the mask from %s ', input_path)
    mask = run_info_parsing.get_mask(input_path)  # read RunInfo.xml
    
    logger.info('Reading SampleSheet from %s', input_path)
    num_lanes, sheet_dict = sample_sheet_parsing.read_sample_sheet(input_path)  # Read SampleSheet.csv
    logger.info(num_lanes)

    sample_project = sample_sheet_parsing.get_sample_project(sheet_dict)

    # Create BCL2FASTQ PBS script
    logger.info('Create BCL2FASTQ pbs script')
    pbs_name = 'BCL_' + project_name + '.pbs'
    logger.info('PBS name: ' + pbs_name)

    create_bcl2fastq_PBS.bcl2fastq_PBS(mask, pbs_name, work_dir, input_path)

    # create BCBIO PBS scripts
    logger.info('Creating BCBIO PBS scripts')
    create_bcbio_PBS.bcbio_loop(sheet_dict, input_path, work_dir, sample_project)
    
    os.chdir(work_dir + '/pbs')
    
    # submit bcl2fastq 
    # create a list with the name of the PBS script
    args = [pbs_name]
    
    # submit the BCL2FASTQ script to batch scheduler
    logger.info('Submitting BCL2FASTQ_PBS')
    BCL2FASTQ_jobid = qsub_dependents.qsub_dependents(args)
    
    # submit the BCBIO scripts once BCL2FASTQ has finished
    logger.info('Submitting BCBIO_PBS')

    list_jobIds = []
    # submit set of BCBIO jobs. A job per sample ID included in the SampleSheet
    for i in range(1, num_lanes):
        script_name = 'runBCBIO_' + str(i) + '.pbs'
        args = [script_name]
        # store the job IDs in a list. They will not get executed until BCL2FASTQ has finished
        list_jobIds.append(qsub_dependents.qsub_dependents(args, BCL2FASTQ_jobid))
         
    logger.info('Done')

