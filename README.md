# Scripts information #
---------------------

*driver.py ->  Main script with calls to the different scripts located at imports directory.

## IMPORTS ##
=======
Directory containing python files with functions used in the main program.
The files included are:

**args.py** -> Parses the arguments of the script. At the moment it only expects the path to a sequencing data.

**checkFinishRun.py** -> Checks when a file has been created indicating the end of a run (Not used at the moment).

**getFastqFiles.py** -> returns a list of all fastq files in a given path (Not used at the moment).a

**xmlparsing.py** -> Parses RunInfo.xml which contains information of the mask required to run BCL2FASTQ

**csvparsing.py** -> Parses a SampleSheet.csv file needed to run BCL2FASTQ. The script expects to find a SampleSheet.csv file within the input directory

**makeProject.py** -> Creates a local project where all the temp files and results will be located before sending them back to the RDF.

**createBCL2FASTQ_PBS.py** -> Creates a PBS script with all the required information to run BCL2FASTQ.

**createBCBIO_PBS.py** -> Create a PBS script with all the required information to run BCBIO. 

**qsub_dependents.py** -> Submits jobs to the queue 


# HOWTO #


```
#!python

driver.py 
```