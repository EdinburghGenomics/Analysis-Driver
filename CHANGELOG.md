Changelog for Analysis-Driver
=============================

0.26 (2019-11-14)
-----------------

- Switch from using clarity REST API to using Reporting LIMS endpoints wherever possible
- Split demultiplexing per lane and allow for non pooling lane to be sequenced in pooling run.
- Fastq-Filterer applied to intermediate fastq files before alignment.
- Project process now use GATK4 for trio check


0.25.3 (2019-09-10)
-------------------

- Fix fastqc command to set the temp directory to lustre.


0.25.2 (2019-08-08)
-------------------

- New toolset for non human variant call with GATK3.8
- Support for Asana API changes


0.25.1 (2019-07-10)
-------------------

- Fixed IDT barcode support
- 'gender' API key renamed to 'sex'


0.25 (2019-06-14)
-----------------

- New pipelines :
  - GATK4 based QC and variants call for human and non-human samples
  - Dragen based variant call for human sample starting from the run processing
- New Features:
  - file program_version.yaml only records version of tools used by the pipeline 
  - Location of genome on the filesystem is provided by the REST API
  - Genome used during Sample processing is uploaded to the REST API
  - Demultiplexing support IDT barcodes
  - Sample Processing will start coverage from run elements is greater than required coverage
  - analysis_driver.log only contains info level logs. New log containing debug level in analysis_driver-debug.log
- Bug fixes:
  - Interop metrics parser doe not uploade NaN
  - Check that all bcl files are present before starting run processing
  - Java tools version can now be determined in cron
  - Increase memory available for trio checking for projects with more than 25 samples 


0.24.1 (2019-03-20)
-------------------

- Bugfix: ensure the bcls expected to exist exists before returning a completed cycle with interop
- Bugfix: make the samplesheet deterministic and only generate it when required


0.24 (2019-03-11)
-----------------

- Alignment during run processing start after the 50 cycles into read2 to speed up alignment metric generation 
- Run picard GC bias detection tool anc calculate new metric to summarise each run element's GC-biais
- Update Picard to version 2.18.23
- Update bcl2fastq to version 2.20
- Bugfix: Initialise analysis_driver_procs at the start to avoid missing embedded entities
- Bugfix: Exceptions raised in Luigi are propagated to the main thread


0.23.1 (2018-11-23)
-------------------

- Updated EGCG-Core to 0.9.1


0.23 (2018-11-09)
-----------------

- Location-independent integration testing with EGCG-Core
- Updated EGCG-Core to v0.9
- Only running trio check on valid projects
- Fixed memory allocation for indel realignment and GenotypeGVCFs
- Eager-loading output config file
- Removed usages of old aggregation
- Overwriting existing fastqs in SampleDataOutput
- Minor fixes: rerunnable stages, pipeline start date, FR insert metrics, PhiX error handling

0.22.1 (2018-06-13)
-------------------

- Compress and index all variants files generated to avoid using GATK own index.


0.22 (2018-06-04)
-----------------

- Upload the source of the processing (run elements for sample, samples for project) to analysis_driver_procs
- Perform variant call for non human in order to get the variant based QC
- Remove Phix reads from data during run processing
- New script for removing Phix reads from already processed data

0.21 (2018-03-15)
-----------------

- Fixed `--resume` option
- Using new Reporting-App aggregation in `ProjectDataset`
- Fixed the toolset_type in `pipelines.demultiplexing`
- New metrics parsed and uploaded:
  - InterOp metrics
  - non-facing read pairs
- Refactored report crawlers to take a single input dir
- Made `relatedness.GenotypeGVCFs` allocate memory dynamically


0.20.1 (2018-02-01)
-------------------

- Add default value for sample that do not have any data
- Fix bug where the sample data threshold was set to be required yield q30 instead of required yield


0.20 (2018-01-30)
-----------------

- Trigger automatic review after sample finishes processing
- Project process duplicated line bugfix
- Update requirements for sample to be ready for processing
- prevent different insert sizes from causing pipeline to crash


0.19.2 (2017-11-27)
-------------------

- Update the field name for required yield/yieldq30/coverage.


0.19.1 (2017-11-22)
-------------------

- Fix temporary directories used by picard.

0.19 (2017-11-16)
-----------------

- In demultiplexing pipeline: align all run elements to their respective default genome. Calculate duplicate rate and coverage
- Each sample processing is recorded as a step in the LIMS

0.18.1 (2017-10-30)
-------------------

- bug fix in run_qc
- buf fix in project process

0.18 (2017-10-11)
-----------------

- Fix bugs in project process
- Rearrange outfile format of GEL relatedness file


0.17 (2017-09-22)
-----------------

- Add Parser for Peddy and Relatedness
- Improvements to integration_test, including more flexible checking of outputs and and retention of data/logs
- Tools versioning: config can now list multiple tool path for multiple version. Toolset config file chose the appropriate version for specific pipeline
- Add md5 files for project process output file
- Original fastq files are now kept when fastq filtering is enabled
- Duplicated unmapped read generated in bcbio are removed
- Small fix to support user-prepared libraries


0.16 (2017-06-22)
-----------------

- Fixed fastq_filterer stats file bug and removed workaround from 0.15.1
- Bug fixes in pipeline stage reporting
- Samples to be passed through variant calling can now be marked `Variant Calling` or `Variant Calling gatk`
- Added genotype relatedness check with [Peddy](https://github.com/brentp/peddy)


0.15.2 (2017-05-29)
-------------------

- Change tools writing to /tmp in Samtools depth and Genotype Validation.


0.15.1 (2017-05-24)
-------------------

- Temporary fix to make sure the fastq_filterer stats file is present to be parsed by RunCrawler.


0.15 (2017-05-18)
-----------------

 - All sample and project processes are now segmented using Luigi
 - Allow filtering/trimming of bad quality runs in demultiplexing 
 - Fix analysis driver --stop and Error handling

0.14.3 (2017-04-28)
-------------------

 - Removed need for SampleSheet to exist for a Run to be picked up


0.14.2 (2017-04-21)
-------------------

 - Fix samplesheet generation when one sample is repeated over multiple lanes
 - Continuously check failed bcl file to make sure they really are failed

0.14.1 (2017-04-19)
-------------------

 - Fix BCL validation bug in previous version

0.14 (2017-04-18)
-----------------

 - Add function to retrieve run metadata from the LIMS in RunDataset, generate the Samplesheet from it
 - Demultiplexing pipeline is now segmented
 - Run processing starts as the first files arrives from the sequencer
 - Bcl validation runs throughout the sequencing
 - Fix bug in species contamination

0.13 (2017-03-07)
-----------------
 - Replace bamtools stats with samtools stats
 - fix bcbio alternative genome version
 - Support for multiple versions of genome per species (configuration file change required)
 - Refactor SampleSheet and RunInfo
 - Support for Seqlab2 sample sheets
 - Archive the output at the end of the process
 - Running pipelines through the [Luigi](http://luigi.readthedocs.io) task runner

0.12 
----
 - BCBio option for variant calling with Freebayes
 - infrastructure for batch processing of entire projects
 - use of new analysis_driver_stages endpoint in Reporting-App
 - split driver into multiple files, contained in pipelines
 - now using egcg_core for notifications
 - updates to Readme.md
 - refactoring throughout

0.11.3
-------
 - Add a try/except statement to catch taxa where taxids are unavailable from the E.T.E TAXDB

0.11.1
-------
- Using [our own fastq filterer](http://github.com/EdinburghGenomics/Fastq-Filterer) after bcl2fastq instead of Sickle
- Added a genome size calculation
- `quiet=True` has been added to Rest API calls in `dataset`, so it should no longer spam the logs in debug mode
- Basic Lims information is now being gathered and pushed to the reporting app early in the pipeline
- Updated EGCG-Core to v0.5.
- Various refactors in contamination_checks, dataset, dataset_scanner, report_crawlers

0.11
-----
- Using our own fastq filterer after bcl2fastq instead of Sickle
- Added a genome size calculation
- quiet=True has been added to Rest API calls in dataset, so it should no longer spam the logs in debug mode
- Basic Lims information is now being gathered and pushed to the reporting app early in the pipeline
- Updated EGCG-Core to v0.5.
- Various refactors in contamination_checks, dataset, dataset_scanner, report_crawlers

0.10.1
--------
- ti/tv and het/hom ration calculation
- bases coverage at 5, 15 and 30X + percentile
