# Scripts information #
---------------------

**driver.py** ->  Main script with calls to the different scripts located at imports directory.

## IMPORTS ##
==============

Directory containing python files with functions used in the main program.
The files included are:

**args.py** -> Parses the arguments of the script. At the moment it only expects the path to a sequencing data.

**checkFinishRun.py** -> Checks when a file has been created indicating the end of a run (Not used at the moment).

**getFastqFiles.py** -> returns a list of all fastq files in a given path (Not used at the moment).a

**xmlparsing.py** -> Parses RunInfo.xml which contains information of the mask required to run `bcl2fastq`.

**csvparsing.py** -> Parses a SampleSheet.csv file needed to run `bcl2fastq`. The script expects to find a SampleSheet.csv file within the input directory

**makeProject.py** -> Creates a local project where all the temp files and results will be located before sending them back to the RDF.

**createBCL2FASTQ_PBS.py** -> Creates a PBS script with all the required information to run `bcl2fastq`.

**createBCBIO_PBS.py** -> Create a PBS script with all the required information to run `bcbio`. 

**qsub_dependents.py** -> Submits jobs to the queue 


# HOWTO #
---------------------

The `driver.py` script expects an path to an input direct:

> python driver.py -i /path/to/sequencing/data

Within that input directory, the scripts will need to find `RunInfo.xml` and `SampleSheet.csv` which contain information related the runs. A working directory, named after the sequencing data, will be created on the local path from 
where the script is executed, therefore it requires writing permissions. Within that new directory, two more directories will be created:

1. Unaligned
2. pbs

The `Unaligned` directory will contain the output of running `bcl2fastq` and `fastqc`, whereas `pbs` will contain the pbs scripts created to run `bcl2fastq` and `bcbio`. In addition, a `bcbio` project
per sample will also be created. 
