# Analysis Driver
---------------------

## Scripts
==========
These can be found in bin/. Currently, these are:

- analysis_driver.py - The main entry point for the pipeline, currently called by ProcTrigger with an input
  data dir.
- run-tests.py - An entry point for pytest unit tests, should the user wish to run the tests outside of
  PyCharm

## Modules
==========

The Analysis Driver consists of several modules, each in turn consisting of several files/functions/classes:

### executor
Classes and functions involved in running Bash commands externally.

- execute() - Takes a list


### reader
Classes that read RunInfo.xml and SampleSheet.csv from the sequencer.

#### run_info
This contains two classes:

- RunInfo - A class that represents an instance of RunInfo.xml. Its purpose is to parse the .xml file and
  construct a Mask object.
- Mask - A class that represents a read mask to be passed to bcl2fastq via self.tostring()

#### sample_sheet
This contains three classes:

- SampleSheet - Reads SampleSheet.csv, validates the barcode lengths contained, and constructs a series of
  SampleProject objects.
- SampleProject - Represents a sample project in SampleSheet.csv, e.g., '10015AT', and contains appropriate
  matching lines from the file. Effectively grouping the sample sheet by sample project. Each SampleProject
  can have its own read length, as long as all its reads are consistent.
- Sample - Represents a line in SampleSheet.csv

### util
Contains various utility classes and functions used by the app, both in __init__ and logger:


- AppLogger - A class that can be mixed in to any class in the Analysis Driver to allow it to write to the
  logging streams set in configuration.
- localexecute - Uses subprocess.Popen to execute arbitrary shell commands
- find_fastqs - Iterates through a directory (intended to be a sample_project folder) and returns a dict of all contained fastq.gz files within, organised per sample.
- setup_bcbio_run - Uses localexecute to run `bcbio -w template` on a given yaml config, a run directory, a
  csv sample file and a list of fastqs

### writer
Contains classes that write various files to be used/executed in the pipeline

#### pbs_writer
This contains four classes that can write PBS-submittable shell scripts:

- PBSWriter - A base class that can open a file for writing, write a PBS header with job parameters
  (`walltime`, `cpus`, `mem`, job name, stdout file and queue id), and save the file.
- BCL2FastqWriter - A PBSWriter that also writes out the bcl2fastq command:
  `bcl2fastq -l INFO --runfolder-dir <run_dir> --output-dir <out_dir> --sample-sheet <sample_sheet> --use-bases-mask <mask>`
- FastqcWriter - Writes the fastqc command:
  `fastqc --nogroup -t 8 -q $(find <input_dir> -name '*.fastq.gz')`
  and also creates an empty `.fastqc_complete` lock file
- BCBioWriter - Exports some Java paths for gatk, and writes out the BCBio command:
  `bcbio <run_yaml> -n <cores> --workdir <workdir>`

#### BCBioCSVWriter
Used to write out a sample file for use with `bcbio -w template`. It takes a reader.SampleSheet object, a
directory to be used in the BCBio run and a directory to search for fastqs. Based on the SampleProject objects
in self.sample_sheet, it writes out a .csv file with the fastq files appropriately grouped.

### config
Defines a class, Configuration, that represents the appropriate environment in the user's .analysisdriver.yaml
config file. Contains configurations for:

- Job submission (currently only PBS)
- Data directories (rdf raw data, input_data, fastqs, BCBio jobs)
- Python interpreter
- Logging (via dictConfig)

### driver
The main 'client' script for the Analysis Driver.

### exceptions
Custom exceptions raised by the pipeline:
- AnalysisDriverError
- ProcessTriggerError


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
