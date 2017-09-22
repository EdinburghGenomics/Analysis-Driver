Changelog for Analysis-Driver
=============================

0.18 (unreleased)
-----------------

- Nothing changed yet.


0.17 (2017-09-22)
-----------------

- Add Parser for Peddy and Relatedness
- Improvements to integration_test, including more flexible checking of outputs and and retention of data/logs
- Tools versioning: config can now list multiple tool path for multiple version. Toolset config file chose the appropriate version for specific pipeline
- Add md5 files for project process output file
- Original fastq files are now kept when fastq filtering is enabled
- Duplicated unmapped read generated in bcbio are removed
- Small fix to support userprepared libraries


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
 - Fix bug in species contamintion

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
