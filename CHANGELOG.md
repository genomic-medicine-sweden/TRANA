# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added

- Added a `generate_master_html` python script that creates `master.html` file containing a table of samples with corresponding pointers to each html output file
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
- Added `Makefile` rules for `precommit`, `lint`, `schema`, `test` as well as a `check` command that runs all checks
- Added `nf-test`
- Added `pre-commit` checks, which includes the `prettier` tool
- Added GitHub Actions workflows for `nf-core` linting checks
- Updated nf-core modules to include nf-tests and added minimal test
- Added emu-combine-output
- Removed authors-field in nextflow.config.
- Added a module to output the result as a phyloseq object

### Fixed

- Fixed unused `NANOPLOT2` to `NANOPLOT_PROCESSED_READS`
- Fixed `Nanoplot` output dirs and prefix in `modules.config`
- Fixed Dockerfile context
- Moved nanostats_unprocessed process execution into seqtype SR if statement
- Conditionally emit nanostats unprocessed/processed to avoid undefined output error when using --seqtype SR
- Fixed a broken configuration file in `modules/local/emu/abundance/meta.yml`
- Fixed failing pulls of some singularity containers by settings singularity.cacheDir

### Changed

- `merge_barcodes_samplesheet.py` can now handle custom barcodes.
- Provided option to `save_merged_reads`
- Provided ability to overwrite files in `publishDir` (`params.publish_dir_overwrite`)
- Cleaned up input preprocessing steps
- Updated `master_template.html` for nanoplot prefix and output changes
- Changed `merge_barcodes_samplesheet.py` container to nf-core (temporary)
- Emit `master.html`
- Changed project name to TACO (`genomic-medicine-sweden/taco`)
- Updated nf-core template to version 3.2.0
- Updated the `ci.yml` config to run test via `nf-test` in addition to the normal test run via cli
- Updated `README.md` with info about running `make install` to gunzip gzipped
  assets and useful `make` commands for developers.

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
