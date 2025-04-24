#!/bin/bash

# author: Olivia Andersson, modified by Katharina Dannenberg
# parameters
# $@: all files to concatenate
# concatenated file will be outpur to stdout, make tsv by collecting to file, e.g. combine_tsv ... > combined.tsv

# input directory
# tsv_dir="/mnt/anvil/Olivia/16s_internutveckling_CG38/gms16s/outdir_internutveckling_gms16s/results"
# tsv_dir=/Users/kadg/Forschung/GMS/16S-Pipeline/gms_16S/results
tsv_dir=$1

# output directory
# output_dir="/mnt/anvil/Olivia/16s_internutveckling_CG38/gms16s/outdir_internutveckling_gms16s"
# output_dir=/Users/kadg/Forschung/GMS/16S-Pipeline/gms_16S/results
# output_dir=$2
# output_dir="/Users/kadg/Forschung/GMS/16S-Pipeline/gms_16S/modules/local/phyloseq"

# Output file name
# output_file="$output_dir/combined-rel-abundance.tsv"

# Find all .tsv files that end with 'rel-abundance.tsv'
# tsv_files=$(find "$tsv_dir" -type f -name '*rel-abundance.tsv')


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

# echo "Combined results saved to $output_file"



