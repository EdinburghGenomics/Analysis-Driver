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
ArrayExecutor that combines multiple Executors into a job array. ClusterExecutor uses `.script_writers` to
write out a PBS, SGE, etc. Bash script and submit it to a queue. All Executor classes have a `join` method
which returns an exit status.

- execute - Takes a list of string Bash commands to execute, and passes it to an appropriate executor. For
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

#### script_writers
This contains classes that can write scripts executable by the shell or by a resource manager:

- ScriptWriter - A base class that can open a file for writing, write commands to it and save.
- PBSWriter - A writer specific to PBSPro. Writes a PBS header to the script upon initialisation (including
  `walltime`, `cpus`, `mem`, job name, stdout file and queue id). Also allows for PBS's job array feature.
- get_script_writer - A function that queries a config object and returns the correct type of writer
  (PBS, SGE, Slurm, etc.) with the correct configurations.


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

#### Notification classes
- Notification - Implements common methods for notifier classes. At the moment, this mainly enforces method
  names and implements formatting of stack traces.
- EmailNotification - Can send HTML-formatted email messages via a configurable SMTP email server. Will
  attempt to reconnect a few times if the connection fails.
- LogNotification - Sends `logging` messages as with other AppLoggers, but also has its own Formatter and
  Handler to log to a central 'notifications' file.


### quality_control
Classes that run further checks on generated output files.

- GenotypeValidation - Implements a workflow using bwa, samtools and gatk to validate called snps against a
  test dataset. Writes a vcf of evaluated genotyping results, and also tries to query the LIMS for an expected
  genotype vcf. If this is successful, it will write a file containing the results of comparing the two vcfs.
- GenderValidation - Quantifies X-chromosome heterozygosity in BCBio's output haplotype vcf. Produces a file
  containing the called gender, to be compared against that suppied in the LIMS.


### reader
Parsers for files supplied and generated during the pipeline.

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


#### mapping_stats_parsers
Functions that pick up QC information from files generated by BCBio, including `bamtools_stats.txt`,
`sort-callable.bed`, `highdepth-stats.yaml` and fastqc reports.

#### demultiplexing_parsers
Functions for picking up QC information from `ConversionStats.xml` and `demultiplexing_stats.xml` (generated
by bcl2fastq) and also files produced by seqtk, containing base counts with adaptors removed.


### util
Contains utility functions:
- bcbio_prepare_samples_cmd - Bash command to execute bcbio_prepare_samples.py for a set of input fastqs,
  merging them. Calls `_write_bcbio_csv`.
- find_files - Convenience function combining `glob.glob` and `os.path.join`. If any files are found for the
  given pattern, returns them.
- find_file - Returns the first result from find_files.
- str_join - Convenience function for calling str.join using `*args`.
- find_fastqs - Uses find_files to return all fastq.gz files. Can limit search to a specific lane.
- find_all_fastqs - Uses `os.walk` to find all `*.fastq.gz` files in a directory.
- write_bcbio_csv - Takes a sample id and a list of fastq files, and writes a csv to be passed to
  bcbio_prepare_samples.

#### bash_commands
Contains functions that build string Bash commands, including:
- bcl2fastq
- fastqc
- seqtk
- bwa mem/samblaster
- bamtools stats
- md5sum
- export of environment variables (for script writers)
- bcbio
- rsync


### report_generation
Contains crawler classes that scan for QC files and pushes them to an external data store:

- RunCrawler
- SampleCrawler


### app_logging
Contains AppLogger, which can be subclassed to implement debug, info, warn, error and critical logging
methods, and get_logger, which generates a logging.logger object and registers it with all active logging
Handlers.


### client
The main 'client' script for the Analysis Driver. Implements argparsing, sets up logging, adds handlers to the
logging configuration and subscribers to the notification centre, calls `dataset_scanner` to find new
datasets, and calls `process_trigger` on the first ready dataset. If the exit status is non-zero or if there
is a stack trace, the dataset will be marked as failed.


### config

#### Configuration 
Singleton, representing the appropriate environment in the user's yaml config file. Contains configurations
for resource manager job submission (currently only PBS), data directories (rdf raw data,
input_data, BCBio jobs), executables (BCBio, fastqc, bcl2fastq, etc.), and logging.

Configuration finds the config file through the environment variable `ANALYSISDRIVERCONFIG`. If this is not
valid, it will then look for `~/.analysisdriver.yaml`. A yaml config must have a `default` environment.
Another environment can be specified with `ANALYSISDRIVERENV`. If another environment is set, Configuration
will use override its `default` parameters with those in the environment. `default`, for example, can  specify
some input data paths, while `testing` instead points to test data. If run with `ANALYSISDRIVERENV=testing`,
the pipeline will override the default and use the testing environment.

Configuration has dict-style [square_bracket] and `self.get()` item retrieval. Can also use `self.query` to
drill down to a parameter.

#### LoggingConfiguration
Contains the pipeline's active logging Handlers. This singleton is set up in `client` using `add_handler`.
`switch_formatter` is used when there is an uncaught Exception to log stacktraces with a blank formatter.
app_logging.get_logger returns a logging.Logger object already registered to the  active Handlers.


### constants
Contains constants for dataset statuses, and database keys/column names used in `rest_communication` and the
report crawlers.


### dataset_scanner
Contains dataset and scanner classes for tracking dataset statuses. Works with the `analysis_driver_procs`
endpoint of the external Rest API

- Dataset - Base class for interacting with the external Rest API.
  - proc_id - Made from the dataset type, name, and start time. Acts as a unique identifier for
    `analysis_driver_procs`.
  - most_recent_proc - Gets the record for the last time the dataset was processed by the pipeline
  - create_process - Adds a new `analysis_driver_proc` for that dataset. Used in `self.change_status`.
  - dataset_status - Queries `self.most_recent_proc` to the get the running status of the dataset. Treats
    `dataset_reprocess` as invisible.
  - change_status - Updates the dataset's most_recent_proc with a new running status. If there is no
    most_recent_proc, then creates one.
  - start - Calls `self.change_status` with 'processing'. Also channges `self.proc_id`, making sure the
    dataset takes a new `analysis_driver_proc`. This way, the Rest API builds up a history of runs for the
    dataset.
  - succeed, fail, abort, reset - Marks a dataset with the relevant status.
  - add_stage - Updates most_recent_proc with a new pipeline stage, including stage name and start time.
  - end_stage - Updates a given stage within most_recent_proc with an end time and exit status.
  - stages - Returns all stages in most_recent_proc that don't have an end date, i.e. currently running.
  - is_ready - `NotImplemented`
- RunDataset - Dataset specific to sequencing runs.
  - rta_complete - Returns presence of RTAComplete.txt
  - is_ready - True if `self.rta_complete` == True, or intermedate_dir is specific in the config.
- SampleDataset - Dataset specific to samples.
  - run_elements
  - read_data - Populates `self.run_elements` with information from the Rest API.
  - runs - Finds all the sequencer runs in `self.run_elements`.
  - amount_data - sums up the q30 bases from `self.run_elements` to get a data volume.
  - data_threshold - Queries the LIMS to get the sample's data threshold, i.e. the minimum volume of data
    required to run variant calling.
  - is_ready - True if `self.amount_data` > `self.data_threshold`.
  - force - Mark the sample for processing even if `self.is ready` is False.

- DatasetScanner - Base class for scanning for datasets.
  - scan_datasets - Uses `self._list_datasets` to list datasets and sort them by status. Returns a dict of
    lists.
  - report - Prints the results of `self.scan_datasets`
  - list_datasets - NotImplemeneted
  - get_dataset - NotImplemented
- RunScanner
  - list_datasets - Lists the contents of the config's run input dir
  - get_dataset - Returns a RunDataset based on a given name and input dir location
- SampleScanner
  - list_datasets - Lists all sample datasets in the Rest API
  - get_dataset - Returns a SampleDataset for the given name



### driver
This contains the main pipeline. `pipeline` takes dataset object and runs demultiplexing, variant calling or
basic QC, depending on the dataset.

- demultiplexing_pipeline - Runs bcl2fastq on sequencer BCL files, runs fastqc, seqtk fqchk and md5sum on
  resulting fastq files and outputs data to the CIFS. Also  outputs the sample sheet and QC files from
  sequencing and demultiplexing.
- variant_calling_pipeline - Finds the relevant fastq files for the (sample) dataset, merges them in the jobs
  folder and runs BCBio variant calling. Also runs fastqc on the merged fastqs, as well as genotype
  validation and gender calling.
- qc_pipeline - Merges fastqs and runs basic bwa alignment. Also runs fastqc and bamtools stats.


### exceptions
Custom exceptions raised by the pipeline:
- AnalysisDriverError


### rest_communication
Contains functions for interacting with the external Rest API. Uses the `requests` module.

- get_documents
- get_document
- post_entry
- put_entry
- patch_entry
- patch_entries - Iterates through payloads and patches each one. Only runs get_documents once.
- post_or_patch - Tries to post an entry and if unsuccessful, prepares the payload and patches instead.

### transfer_data
Functions for transferring data between the CIFS and local temporary storage.

- prepare_run_data - Runs `_transfer_run_to_int_dir` if an intermediate processing directory is set.
- prepare_sample_data - Finds all fastq files for all run elements in the dataset object.
- find_fastqs_for_sample - First looks in `jobs` for fastqs. If none there, looks in `input_dir`. If
  `input_dir` is remote, runs an rsync to `jobs` and checks again.
- transfer_run_to_int_dir - rsyncs the data across until it is
  known that the sequencing is complete.
- create_links_from_bcbio - Uses etc/output_files.yaml to find output files to be picked up and symbolically
  links them to a 'to be outputted' directory in the jobs folder.
- output_run_data - rsyncs the data across.
- output_sample_data - rsyncs the data across, but also queries the LIMS for the user's supplied sample ID.



# HOWTO #
---------------------
The Analysis Driver is run from bin/edingen_analysis_driver.py:

    python bin/edingen_analysis_driver.py <args>

Valid arguments include:
- --debug - Run with log levels set to `DEBUG`
- --report - Run `scan_datasets` only. Will not list (potentially numerous) completed/failed datasets.
- --report-all - As above, but report all datasets.
- --skip <datasets...> - Mark the specified datasets as `complete`.
- --reset <datasets...> - Mark the specified datasets with `reprocess`.
- --abort <datasets...> - Mark specified datasets as `aborted`. These will not be picked up by the dataset
  scanner

If no args are given, one new dataset will be processed. Run recurrently through Cron to periodically scan
datasets and kick off one

To ignore datasets (or indeed any non-dataset directories) in the configured location, a `.triggerignore` file
can be written in the same place, where each line is the folder name. Any directories listed in this file will
not be picked up by the scanner.

Within any input dataset, there will need to be a `RunInfo.xml` and a `SampleSheet.csv`. These contain
information about the sequencing runs. A working directory will be created in the (configured) jobs folder
where intermediate and output data is kept.
