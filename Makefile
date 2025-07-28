# ==============================================================================
# Brief explanation of Makefile syntax
# ==============================================================================
# Since Makefiles are a bit special, and there are not a lot of easy to read
# tutorials out there, here follows a very brief introduction to the syntax in
# Makefiles.
#
# Basics
# ------------------------------------------------------------------------------
#
# Makefiles are in many ways similar to bash-scripts, but unlike bash-scripts
# they are not written in a linear sequential fashion, but rather divided into
# so called rules, that are typically tightly connected to specific output file
# paths or file path pattern.
#
# Firstly, Makefiles should be named `Makefile` and be put inside a folder
# where one wants to run the make command.
#
# The rule syntax
# ------------------------------------------------------------------------------
#
# The syntax of a rule is at its core very simple:
#
# output file(s)/target(s) : any dependent input files
#     commands to produce the output files from input files
#     possibly more commands
#
# So, the significant part of the rule syntax is the ':'-character, as well as
# the indentation of the rows following the head of the rule, to indicate the
# "recipe" or the commands to produce the outputs of the rule.
#
# A rule can also have just a name, and no concrete output file. That is, it
# would have the form:
#
# name_of_the_rule : any dependent input files
#     one commands
#     more commands
#
# Now, there is one big caveat here, related to our scripts: Make will rebuild
# a target as soon as any of its inputs are updated, or have a newer timestamp
# than the target. This is typically not desired for us, since we might have
# files unpacked with zip with all kinds of different dates, and for a one-run
# installation, we are mostly interested in wheter an output file already
# exists or not, not so much about timestamps.
#
# To change so that Make only cares about whether files exist, and not timestamps,
# one can add a | character before those input files, like so:
#
# name_of_the_rule : input files where timestamp matters | input files where only existence matters
#     one commands
#     more commands
#
# Of course, one can have everything on the right side of the | character, so:
#
# name_of_the_rule : | input files where only existence matters
#     one commands
#     more commands
#
# Running Makefiles
# ------------------------------------------------------------------------------
#
# To run a named rule, you would then do:
#
# $ make name_of_the_rule
#
# For normal rules, you would hit:
#
# $ make <outputfile>
#
# Tip: Type "make" and hit tab twice in the shell, to show available targets to
# run in the current makefile.
#
# Special variables and more
# ------------------------------------------------------------------------------
#
# Inside the commands of the rules, one can write pretty much bash, especially
# if setting `SHELL := /bin/bash` in the beginning of the Makefile which we have
# done below.
#
# There are a some differences though:
#
# 1. Variable syntax using only a single $-character refer to MAKE-variables.
#    Thus, to use variables set in bash, you have to use double-$.
#    So:
#    echo $(this-is-a-make-variable) and $$(this-is-a-bash-variable)
#
#    This also goes for command substitution, so rather than:
#
#    	echo "Lines in file:" $(wc -l somefile.txt)
#
#    ... you would write:
#
#    	echo "Lines in file:" $$(wc -l somefile.txt)
#
# 2. Makefiles also use some special variables with special meaning, to access
#    things like the output and input files:
#
#    $@   The output target (In case of multiple targets, the same command will
#         be run once for every output file, feeding this variable with only a
#         single one at a time)
#
#    $<   The first input file or dependency
#
#    $^   ALL input files / dependencies.
#
# 3. Another speciality is that adding a '@' character before any command, will
#    stop this command from being printed out when it is executed (only its output
#    will be printed).
#
# That's all for now. For the official docs, see:
# https://www.gnu.org/software/make/manual/make.html
#
# Some more caveats
# ------------------------------------------------------------------------------
#  Here are some things that might not be obvious at first:
#
#  To create a rule with a single output but many dependencies, you can add
#  these dependencies on multiple lines by using the continuation character \
#  like so:
#
#  name_of_rule: dependency1 \
#  		dependency2 \
#  		dependency3 \
#  		dependency4 \

# ==============================================================================
# Various definitions
# ==============================================================================
# Make sure bash is used as shell, for consistency and to make some more
# advanced scripting possible than with /bin/sh
SHELL := /bin/bash

# Define path variables
SCRIPT_DIR := $(shell pwd)
ASSETS_DIR := $(shell realpath $(SCRIPT_DIR)/assets/)
CONTAINER_DIR := $(realpath $(SCRIPT_DIR)/container/)
# The root folder where the pipeline is currently located. To be mounted into
# the Singularity containers below.
MNT_ROOT := /$(shell readlink -f . | cut -d"/" -f2)
INSTALL_LOG := "$(SCRIPT_DIR)/.install.log"

define log_message
	@echo "--------------------------------------------------------------------------------" | tee -a $(INSTALL_LOG);
	@echo "$$(date "+%Y-%m-%d %H:%M:%S"): $1" | tee -a $(INSTALL_LOG);
	@echo "--------------------------------------------------------------------------------" | tee -a $(INSTALL_LOG);
endef

# ==============================================================================
# Main rules
# ==============================================================================

install: assets/databases/emu_database/species_taxid.fasta assets/databases/emu_database/taxonomy.tsv assets/databases/krona/taxonomy/taxonomy.tab

assets/databases/emu_database/species_taxid.fasta: assets/databases/emu_database/species_taxid.fasta.gz
	zcat $< > $@

assets/databases/emu_database/taxonomy.tsv: assets/databases/emu_database/taxonomy.tsv.gz
	zcat $< > $@

assets/databases/krona/taxonomy/taxonomy.tab: assets/databases/krona/taxonomy/taxonomy.tab.gz
	zcat $< > $@

schema:
	nf-core pipelines schema build

precommit:
	pre-commit run --all-files

nf-core-lint:
	nf-core pipelines lint

nextflow-lint:
	nextflow lint -o concise main.nf

lint: nextflow-lint nf-core-lint

check: precommit lint

nf-test:
	nf-test test

test-cli-fastq:
	nextflow  \
		-log $$(pwd)/nextflow.log \
		run main.nf \
		-profile singularity,test \
		--outdir results \
		--db $$(pwd)/assets/databases/emu_database \
		--merge_fastq_pass $$(pwd)/assets/test_assets/ci

test-cli-samplesheet:
	nextflow  \
		-log $$(pwd)/nextflow.log \
		run main.nf \
		-profile singularity,test \
		--outdir results \
		--db $$(pwd)/assets/databases/emu_database \
		--input https://raw.githubusercontent.com/genomic-medicine-sweden/test-datasets/refs/heads/16s/samplesheet.csv

test: nf-test test-cli-fastq test-cli-samplesheet

check-and-test: check test
