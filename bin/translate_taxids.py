#!/usr/bin/env python

import argparse
import pandas as pd
import sys

VERSION = "v0.0.1"

def get_best_tax_label(row):
    """Return the best available taxonomic label from left to right."""
    for level in ['species', 'genus', 'family', 'order', 'class', 'phylum', 'clade', 'superkingdom', 'subspecies', 'species subgroup', 'species group']:
        val = row.get(level, '')
        if pd.notna(val) and str(val).strip() != '':
            return val.strip()
    return 'Unknown'

def main(tsv_path, taxonomy_path, output_path):
    df = pd.read_csv(tsv_path, sep='\t', header=0)
    original_headers = df.columns.tolist()

    taxonomy_df = pd.read_csv(taxonomy_path, sep='\t', dtype=str).set_index('tax_id')
    tax_id_to_label = {
        tax_id: get_best_tax_label(row)
        for tax_id, row in taxonomy_df.iterrows()
    }

    new_headers = [
        tax_id_to_label.get(col.strip(), col) if col.strip().isdigit() else col
        for col in original_headers
    ]
    df.columns = new_headers

    df.to_csv(output_path, sep='\t', index=False)
    print(f"Output written to {output_path}")

if __name__ == '__main__':
    # First, only parse --version if it’s present
    parser = argparse.ArgumentParser(add_help=False)
    parser.add_argument('--version', action='store_true')
    args, remaining_args = parser.parse_known_args()

    if args.version:
        print(f"translate_taxids.py version {VERSION}")
        sys.exit(0)

    # Now parse the full arguments including required positional ones
    parser = argparse.ArgumentParser(description="Translate tax_id columns to taxonomic names.")
    parser.add_argument("tsv_file", help="Input TSV file with tax_id headers")
    parser.add_argument("taxonomy_file", help="Taxonomy file with tax_id and taxonomic ranks")
    parser.add_argument("output_file", help="Output TSV file with translated headers")
    args = parser.parse_args(remaining_args)

    main(args.tsv_file, args.taxonomy_file, args.output_file)
