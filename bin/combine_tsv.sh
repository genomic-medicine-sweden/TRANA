#!/bin/bash

# author: Olivia Andersson, modified by Katharina Dannenberg
# version: 1.0
# parameters:
# $@: all files to concatenate
# concatenated file will be output to stdout, make tsv by collecting to file, e.g. combine_tsv ... > combined.tsv

# Initialize a flag to add the header only once
header_added=false

for tsv in $@; do
# echo "$tsv"
    if [ "$header_added" = false ]; then
        head -n 1 "$tsv" | awk '{print $0 "\tSource_File"}'
        tail -n +2 "$tsv" | awk -v fname=$(basename "$tsv") '{print $0 "\t" fname}'
        header_added=true
    else
        tail -n +2 "$tsv" | awk -v fname=$(basename "$tsv") '{print $0 "\t" fname}'
    fi
done


