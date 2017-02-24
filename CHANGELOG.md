Changelog for Analysis-Driver
=============================

0.13 (Unreleased)
-----------------
 - Replace bamtools stats with samtools stats

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
Using our own fastq filterer after bcl2fastq instead of Sickle
Added a genome size calculation
quiet=True has been added to Rest API calls in dataset, so it should no longer spam the logs in debug mode
Basic Lims information is now being gathered and pushed to the reporting app early in the pipeline
Updated EGCG-Core to v0.5.
Various refactors in contamination_checks, dataset, dataset_scanner, report_crawlers

0.10.1
--------
- ti/tv and het/hom ration calculation
- bases coverage at 5, 15 and 30X + percentile
