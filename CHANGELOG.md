# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [v0.4.2]

### Added

### Fixed

- Fixed `filtlong` bug that crashes pipeline when all reads are filtered out with `.branch` operator
- Fixed fastqc data missing in multiqc report
- No v in version in nextflow.config

### Changed

## [v0.4.1]

### Added

### Fixed

- Fixed version in nextflow.config

### Changed

## [v0.4.0]

### Added

## [v0.4.0]

### Added

- Added a `generate_master_html` python script that creates `master.html` file
  containing a table of samples with corresponding pointers to each html output
  file
- Added repective `GENERATE_MASTER_HTML` process
- Added `cmd.config`
- Added `params.trace_timestamp` to `nextflow.config`
- Added `changelog_update_reminder` GA workflow
- Added optional ability to save merged reads
- Added capability to emit nanoplot nanostats text file
- Added clear commenting for modules and subworkflows
- Added stubs for all local modules
- Added `when` operator for nanoplot process
- Added workflow that publishes docker image for python-related processes
- Added contributors
- Added versioning to `merge_barcodes_samplesheet.py`
- Added output argument to `generate_master_html.py`
- Added prefix to `master.html` filename
- Added `seqtk_sample` process
- Added `Makefile` rules for `precommit`, `lint`, `schema`, `test`,
  `test-cli-fastq`, `test-cli-samplesheet`, as well as a `check` command that
  runs all checks
- Added `nf-test`
- Added `pre-commit` checks, which includes the `prettier` tool
- Added GitHub Actions workflows for `nf-core` linting checks
- Updated nf-core modules to include nf-tests and added minimal test
- Added emu-combine-output
- Removed authors-field in nextflow.config.
- Added a module to output the result as a phyloseq object
- Added figure for pipeline
- Added an nf-test for short-read, paired-end data
- Added tests for the combinations of adapter trimming and quality filtering parameters.
- Added test for downsampling
- Added running of `nextflow lint` to CI
- Added module for heatmap for emu likelihood file
- Added barplots for control comparison
- Added module for translation of tax-ids into scientific name for heatmap generation
- Logs for local modules
- result directory for seqtk
- Added a config that can be used for SGE

### Fixed

- Readded `NANOPLOT2` as `NANOPLOT_PROCESSED_READS`
- Fixed Dockerfile context
- Moved nanostats_unprocessed process execution into seqtype SR if statement
- Conditionally emit nanostats unprocessed/processed to avoid undefined output error when using --seqtype SR
- Fixed a broken configuration file in `modules/local/emu/abundance/meta.yml`
- Fixed failing pulls of some singularity containers by settings singularity.cacheDir
- Fixed broken nanoplot code. Readded prefixes to the module.
- Fixed broken `--help` screen by adding missing nf-schema plugin and validation config.
- Fixed broken e-mail sending, where a (corrupt) e-mail would be sent at the start of a pipeline rather than at the end.
- Fixed a pipeline ordering error in the main workflow where Nanoplot outputs were accessed before executing the module.
- Fixed various broken process calls in the main Trana workflow found when adding tests for parameter combinations.
- Updated schema
- Fixed input channel for seqtk
- Fixed names for containers
- Fixed seqtk links

### Changed

- `merge_barcodes_samplesheet.py` can now handle custom barcodes.
- Provided option to `save_merged_reads`
- Provided ability to overwrite files in `publishDir` (`params.publish_dir_overwrite`)
- Cleaned up input preprocessing steps
- Updated `master_template.html` for nanoplot prefix and output changes
- Changed `merge_barcodes_samplesheet.py` container to nf-core (temporary)
- Emit `master.html`
- Changed project name to TRANA (`genomic-medicine-sweden/trana`)
- Updated nf-core template to version 3.2.0
- Updated the `ci.yml` config to run test via `nf-test` in addition to the normal test run via cli
- Updated `README.md` with info about running `make install` to gunzip gzipped
  assets and useful `make` commands for developers.
- Updated workflow image
- Made phyloseq an option to protect against errors when using other databases than emu´s premade.
- Added peoples github names to README for credit
- Added nanoplot unprocessed reads and processed reads to multiqc. Short info about this in README.
- Added steps in `PIPELINE_INITIALISATION` from nf-core template, and replaced workflow.onComplete code
  with call to `PIPELINE_COMPLETION` subworkflow
- Updated images for github
- Added versions for local scripts
- Changed `PIPELINE_INITIALISATION` subworkflow to be modularised
- Changed `bin/merge_barcodes_samplesheet.py` to take csv as input instead of tsv as stated in `README` (safer option)
- Name of pipeline was change to TRANA

## [v0.1.0]

### Added

- Support for Nanopore reads. Reads can be in barcode directories or premerged into samplespecific fastq-files.
- Krona
- emu
- FiltLong
- Porchop_abi
- NanoPlot
- Multiqc
- fastqc
- Adaptation to nextflow
- All tools are containerised with singularity
