# Analysis Driver
[![travis](https://img.shields.io/travis/EdinburghGenomics/Analysis-Driver/master.svg)](https://travis-ci.org/EdinburghGenomics/Analysis-Driver)
[![landscape](https://landscape.io/github/EdinburghGenomics/Analysis-Driver/master/landscape.svg)](https://landscape.io/github/EdinburghGenomics/Analysis-Driver)
[![GitHub issues](https://img.shields.io/github/issues/EdinburghGenomics/Analysis-Driver.svg)](https://github.com/EdinburghGenomics/Analysis-Driver/issues)

Dependencies:
- [Reporting-App](https://github.com/EdinburghGenomics/Reporting-App) - Rest API for storing QC data and
  tracking the running of pipeline stages
- [EGCG-Core](https://github.com/EdinburghGenomics/EGCG-Core) - Interaction with the Rest API and Clarity
  Lims, HSM data archiving, notifications, and resource manager job execution

## Scripts
These can be found in `bin/`. Currently, these are:

- edingen_analysis_driver.py - The main entry point for running the pipeline
- run_qc.py - Entry point for rerunning QC data generation for a given dataset
- send_data.py - Entry point for rerunning the QC crawlers for a given dataset, pushing/re-pushing QC data
  to the Rest API


## Architecture
The Analysis Driver consists of several modules, each in turn consisting of several files/functions/classes:

### quality_control
Classes that run checks on output files generated from the main pipeline.

- BCLValidator - Runs `gzip -t` on BCL files in a directory and records exit statuses in a csv record
- Relatedness - Checks relatedness for multiple samples using `vcftools relatedness2`
- BadCycleTileDetector - Inspects a run's InterOp files for low-quality tiles and cycles to pass to
  [fastq_filterer](https://github.com/EdinburghGenomics/Fastq-Filterer)
- GenotypeValidation - Uses bwa, samtools and gatk to validate called snps against a test dataset. Compares
  a sample's genotype with an expected, queries the LIMS for an expected genotype vcf, and writes a file containing the
  results of comparing the observed and expected vcfs
- GenderValidation - Quantifies X-chromosome heterozygosity in BCBio's output haplotype vcf. Produces a file
  containing the called gender, to be compared against the gender suppied in the Lims
- FastqScreen - Checks fastqs for sample contamination using [fastqscreen](http://www.bioinformatics.babraham.ac.uk/projects/download.html#fastqscreen)
- Blast - Checks a fastq file for contamination using NCBI Blast
- VerifyBamID - Checks a Bam file for species contamination using [VerifyBamID](http://genome.sph.umich.edu/wiki/VerifyBamID)
- VCFStats - Checks for contamination using VCF metrics from [rtg vcfstats](https://github.com/RealTimeGenomics/rtg-tools)
- SamtoolsDepth - Uses samtools to calculate median coverage in a Bam file
- WellDuplicates - Runs [well_duplicates](https://github.com/EdinburghGenomics/well_duplicates), a tool for
  predicting duplicate reads in the flowcell

### reader
Parsers for files supplied and generated during the pipeline.

- run_info
  - RunInfo - Represents an instance of RunInfo.xml. Reads in the xml and constructs a Mask object from the
    `<Reads>` element.
  - Mask - Represents a read mask to be passed to bcl2fastq via sample_sheet.generate_mask.
- mapping_stats_parsers - Functions that pick up QC information from files generated by BCBio, including
  `bamtools_stats.txt`, `sort-callable.bed`, `highdepth-stats.yaml` and fastqc reports.
- demultiplexing_parsers - Functions for picking up QC information from `ConversionStats.xml` and
  `demultiplexing_stats.xml` (generated by bcl2fastq) and also files produced by seqtk, containing base counts
  with adaptors removed.

### util
Miscellaneous utilies and a collection of functions that build string Bash commands for running bcl2fastq,
fastqc, BCBio, etc.

### report_generation
Contains crawlers that scan for QC files, parses them using functions from `reader` and pushes data to the
Rest API.

### client
The main 'client' script for the Analysis Driver. Sets up argparsing, logging and notifications, calls
`dataset_scanner` to find new datasets, runs the first ready dataset, and marks it as successful or failed
depending on the exit status.

### config
Contains several config objects, which read various files in search paths and in `etc/`.

### dataset
When a pipeline is running, data about the pipeline run is handled by a `Dataset` object. The Dataset
interacts with external data sources and controls notifications. The dataset in turn contains a
`MostRecentProc` object, which handles the stages of the pipeline run.

`Dataset.dataset_status` gets the running status of the dataset - used by `dataset_scanner` to determine
whether to process the dataset. `dataset_reprocess` is treated as invisible, and will show up as either
`dataset_new` or `dataset_ready` depending on the dataset's ready status.

The dataset is able to stop an already-running pipeline by sending it a SIGUSR1 or SIGUSR2. SIGUSR2 will stop
a pipeline where it is, while SIGUSR1 will cause the Luigig task runner to stop scheduling any new taks,
stopping the pipeline cleanly.

RunDataset is able to query the Lims and the run's `RunInfo.xml` for data. It can also write this data to a
sample sheet to be passed to bcl2fastq.

SampleDataset can query the Rest API for run elements and the Lims for the sample's Yield for Quoted Coverage.

ProjectDataset can retrieve the species and genome version of a project.

### notification
Subclasses the EGCG-Core `NotificationCentre`, adding pipeline-specific messaging for `start_pipeline`,
`end_pipeline`, `start_stage`, `end_stage` and `crash_report`.

### dataset_scanner
Dataset scanners scan the Rest API, returning datasets sorted by dataset status. This is then used by the
`client` to run a new dataset, and can also be printed to stdout for reporting.

### pipelines
This contains the various logic workflows that Analysis-Driver is capable of.

- pipeline - Takes a dataset object, decides which pipeline to use and runs it.
- demultiplexing - Runs bcl2fastq, runs QC and filtering on the resulting data, and outputs it. Also outputs
  extra files produced during sequencing and demultiplexing.
- var_calling - Finds the relevant fastq files for the SampleDataset, merges them in the job folder
  and runs bam_file_production and GATK variant calling.
- bam_file_production - Merges fastqs and runs bwa alignment. Also runs fastqc and samtools stats.
- qc - Runs bam_file_production and data output.
- bcbio - Merges fastqs and runs variant calling via BCBio. Can use GATK or Freebayes.

All pipelines are segmented into [Luigi](http://luigi.readthedocs.io) tasks, allowing the pipeline to restart
a failed dataset from the last successful stage. The classes in `quality_control` are also implemented as
Luigi tasks.

### transfer_data
Functions for finding, outputting and archiving output files.

## Usage
The Analysis Driver is run from bin/edingen_analysis_driver.py:

    python bin/edingen_analysis_driver.py (--run | --sample | --project) [--debug] <action>

Valid actions are:
- --report - Scan for datasets only. Will not list (potentially numerous) completed/failed datasets
- --report-all - As above, but report all datasets
- --skip <datasets...> - Mark a dataset as `complete`
- --abort <datasets...> - Mark a dataset as `aborted` - these will not be picked up by the dataset scanner
- --reset <datasets...> - Mark a dataset for re-running from the start
- --resume <datasets...> - Mark a dataset for re-running from the last completed stage
- --force <datasets...> - Mark a sample for processing, even if below the data threshold
- --stop <datasets...> - Stop the PID of a running dataset
- --soft-stop <datasets...> - Stop scheduling new tasks for a running dataset

If no action is given, one new dataset will be processed. Run recurrently through Cron to periodically
scan datasets and kick off one at a time.

To ignore datasets (or indeed any non-dataset directories) in the configured location, a `.triggerignore` file
can be written in the same place, where each line is the folder name. Any directories listed in this file will
not be picked up by the scanner.

## Integration tests
There is a set of end-to-end pipeline tests to be run on compute infrastructure. These require test input
data, a config file describing expected output data, and a Docker image of
[Reporting-App](https://github.com/EdinburghGenomics/Reporting-App), installed under the name
`egcg_reporting_app`.
