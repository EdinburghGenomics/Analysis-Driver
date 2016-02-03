# Analysis Driver
---------------------

## Scripts
==========
These can be found in bin/. Currently, these are:

- edingen_analysis_driver.py - The main entry point for the pipeline, currently called by ProcTrigger with an input
  data dir.
- run-tests.py - An entry point for pytest unit tests, should the user wish to run the tests outside of
  PyCharm


## Modules
==========

The Analysis Driver consists of several modules, each in turn consisting of several files/functions/classes:

### executor
Classes and functions involved in running Bash commands externally. Executor classes that run locally include
a non-threaded Executor, a threaded StreamExecutor that streams its stdout/stderr as it runs, and an
ArrayExecutor that combines multiple Executors into a job array. ClusterExecutor calls
`analysis_driver.writer` to write out a PBS, SGE, etc. Bash script and submit it to a queue. All Executor
classes have a `join` method which returns an exit status.

- `execute()` - Takes a list of string Bash commands to execute, and passes it to an appropriate executor. For
  a single command, pass in a list of one string. By default, it will use ClusterExecutor unless
  `env` is 'local', in which case a local executor will be used.
- ClusterExecutor - Converts a Bash command into a `subprocess.Popen` object, which is executed when
  `self.join` is called.
- StreamExecutor - Subclasses Executor and `threading.Thread`, making its Popen run in a sub-thread. To
  use, generate the object then either call `self.run` or call `self.start` and then ideally `self.join`
  later.
- ArrayExecutor - Takes a list of Bash commands and generates an Executor object for each. In `self.run`, it
  iterates over these and calls `start` and `join` for each one.
- ClusterExecutor - Takes a list of Bash commands. Upon creation, it calls `writer.get_script_writer` and
  writes a Bash script. `prelim_cmds` may be specified here to, e.g, export Java paths prior to commencing the
  job array. `self.start` and `self.join` executes a `qsub` Bash command for its script.

### notification
Classes able to send messages and notifications over various media. To use, import
`notification_center.default` and use its subscriber methods: `start_pipeline`, `end_pipeline`, `start_stage`
and `end_stage`.

#### NotificationCenter
Singleton for import by other modules. The `client` module calls `self.add_subscribers` according to its
configuration, which appends a list of subscribers. Also implements `self._pass_to_subscriber`, which passes
on any called methods to its subscribers. Therefore, calling `notification_center.some_method` will try to
call `some_method` on each of its subscribers. If the subscriber does not implement the method, a `debug`
message will be logged.

#### Notification
Implements common methods for notifier classes. At the moment, this mainly enforces method names and
implements formatting off stack traces.

#### EmailNotification
Has methods for sending messages via a configurable SMTP email server. The message is formatted by a Jinja2
template and sent by `smtplib`. If there is an error connecting to the server, it will attempt to reconnect
up to 3 times.

#### LogNotification
Sends messages via `logging` as with other AppLoggers, but also has its own Formatter and Handler to log to a
central notifications file.

### quality_control

#### GenoTypeValidation
This implements a workflow using bwa, samtools and gatk to validate called snps against a test dataset.


### reader
Classes that read RunInfo.xml and SampleSheet.csv from the sequencer.

#### run_info
This contains two classes:

- RunInfo - A class that represents an instance of RunInfo.xml. Its purpose is to parse the .xml file and
  construct a Mask object.
- Mask - A class that represents a read mask to be passed to bcl2fastq via self.tostring()

#### sample_sheet
This contains three classes:

- SampleSheet - Reads SampleSheet_analysis_driver.csv, validates the barcode lengths contained, and constructs
a series of SampleProject objects.
- SampleProject - Represents a sample project in SampleSheet.csv, e.g., '10015AT', and contains appropriate
  matching lines from the file. Effectively grouping the sample sheet by sample project. Each SampleProject
  can have its own read length, as long as all its reads are consistent.
- Sample - Represents a line in SampleSheet.csv

Also present is transform_sample_sheet, which takes SampleSheet.csv, transforms it to a format valid for
bcl2fastq, and writes it out as SampleSheet_analysis_driver.csv. It is this transformed sample sheet that is
passed to SampleSheet.


### util
Contains utility functions:
- bcbio_prepare_samples - Executes bcbio_prepare_samples.py for a seet of input fastqs, merging them.
- setup_bcbio_run - Runs `bcbio -w template` on a given yaml config, run directory, csv sample file and a list
  of fastqs.
- transfer_output_files - Uses `shutil.copyfile` to transfer pipeline output files to a configured directory.
  Also renames output files according to configuration (see etc/example_analysisdriver.yaml).


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


### app_logging
Contains AppLogger, which can be subclassed to implement debug, info, warn, error and critical logging
methods, and get_logger, which generates a logging.logger object and registers it with all active logging
Handlers.


### client
The main 'client' script for the Analysis Driver. Implements argparsing, sets up logging, adds handlers to the
logging configuration and subscribers to the notification centre, calls `dataset_scanner` to find new
datasets, and calls `process_trigger` on the first ready dataset. If the exit status is non-zero or if there
is a stack trace, a 'failed' lock file will be set for the dataset.


### config

#### Configuration 
Singleton, which represents the appropriate environment in the user's yaml config file. Contains
configurations for resource manager job submission (currently only PBS), data directories (rdf raw data,
input_data, BCBio jobs), executables (BCBio, fastqc, bcl2fastq, etc.), renaming of output files and sample
sheet fields, and logging.

Configuration finds the config file through the environment variable `ANALYSISDRIVERCONFIG`. If this is not
valid, it will then look for `~/.analysisdriver.yaml`. A yaml config must have a `default` environment. Upon
creation of the singleton, it will take this default and, if necessary, override its parameters with those
from another environment specified in `ANALYSISDRIVERENV`. `default`, for example, can  specify some input
data paths, while `testing` instead point to test data. If run with `ANALYSISDRIVERENV=testing`, the pipeline
will override the default and use testing environment.

Configuration has dict-style [square_bracket] and `self.get()` item retrieval. Can use `self.query` to drill
down to a parameter. Can also return a string representation of its content for reporting via `self.report`.

#### LoggingConfiguration
Contains the pipeline's active logging Handlers. This singleton is set up in `client` using `add_handler`.
`switch_formatter` is used when there is an uncaught Exception to log stacktraces with a blank formatter.
app_logging.get_logger returns a logging.Logger object already registered to the  active Handlers.


### dataset_scanner
Provides functions for identifying datasets in the input directory. The scanner uses lock files alongside
(usually) the datasets. This may become a database in the future.

- dataset_status - Locates a dataset's lock file and returns a 'status', based on its extension. If there is
  no lock file, the dataset is treated as 'new'. Only one lock file may exist for a dataset at once.
- scan_datasets - Uses dataset_status to report an all datasets in the input directory.
- reset - Removed any lock file for a dataset, effectively making it 'new' again.
- switch_status - `reset`s a dataset, then touches a lock file for the desired status. Used for skipping or
  ignoring datasets.


### driver
This contains the main pipeline. `pipeline` takes a path to an input dataset, and runs bcl2fastq, fastqc,
genotype validation and variant calling with BCBio. If no sample sheet is found, or if RunInfo.xml has no
barcode reads, an alternate phiX pipeline will be used, only running bcl2fastq and fastqc. Returns an exit
status.


### exceptions
Custom exceptions raised by the pipeline:
- AnalysisDriverError
- ProcessTriggerError


### process_trigger
- trigger - If an intermediate processing directory is set, rsync the data across until it is known that the
  sequencing is complete. Then run `driver.pipeline` and return an exit status.



# HOWTO #
---------------------
The Analysis Driver is run from bin/edingen_analysis_driver.py:

    python bin/edingen_analysis_driver.py <args>

Valid arguments include:
- --debug - Run with log levels set to `DEBUG`
- --report - Run `scan_datasets` only. Will not list (potentially numerous) completed/failed datasets.
- --report-all - As above, but report all datasets.
- --skip <datasets...> - Set a `complete` lock file for the specified datasets.
- --reset <datasets...> - Reset lock files for specified datasets.
- --abort <datasets...> - Mark specified datasets as aborted. These will not be picked up by the dataset
  scanner

If no args are given, one new dataset will be processed. Run recurrently through Cron to periodically scan
datasets and kick off one

To ignore datasets (or indeed any non-dataset directories) in the configured location, a `.triggerignore` file
can be written in the same place, where each line is the folder name. Any directories listed in this file will
not be picked up by the scanner.

Within any input dataset, there will need to be a `RunInfo.xml` and a `SampleSheet.csv`. These contain
information about the sequencing runs. If no sample sheet is found or RunInfo.xml has no barcode reads, a
smaller phiX pipeline will be used. A working directory will be created in the (configured) jobs folder where
intermediate and output data is kept.
