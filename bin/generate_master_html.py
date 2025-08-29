#!/usr/bin/env python

"""Generate a master html template."""

import os
import re
import argparse
import pandas as pd
from jinja2 import Template
from datetime import datetime

description = '''
------------------------
Title: generate_master_html.py
Date: 2024-12-16
Author(s): Ryan Kennedy
------------------------
Description:
    This script creates master html file that points to all html files that were outputted from EMU.

List of functions:
    get_sample_ids, generate_master_html.

List of standard modules:
    csv, os, argparse.

List of "non standard" modules:
    pandas, jinja2.

Procedure:
    1. Get sample IDs by parsing samplesheet csv.
    2. Render html using template.
    3. Write out master.html file.

-----------------------------------------------------------------------------------------------------------
'''

usage = '''
-----------------------------------------------------------------------------------------------------------
Generates master html file that points to all html files.
Executed using: python3 ./generate_master_html.py -i <Input_Directory> -o <Output_Filepath>
-----------------------------------------------------------------------------------------------------------
'''

parser = argparse.ArgumentParser(
                description=description,
                formatter_class=argparse.RawDescriptionHelpFormatter,
                epilog=usage
                )
parser.add_argument(
    '-v', '--version',
    action='version',
    version='%(prog)s 0.2.0'
    )
parser.add_argument(
    '-c', '--csv',
    help='input samplesheet csv filepath',
    metavar='SAMPLESHEET_CSV_FILEPATH',
    dest='csv',
    required=True
    )
parser.add_argument(
    '-m', '--html',
    help='input master html template filepath',
    metavar='MASTER_HTML_TEMPLATE_FILEPATH',
    dest='html',
    required=True
    )
parser.add_argument(
    '-t', '--timestamp',
    help='pipeline execution timestamp',
    metavar='PIPELINE_EXECUTION_TIMESTAMP',
    dest='timestamp',
    required=True
    )
parser.add_argument(
    '-o', '--output',
    help='output filepath',
    metavar='OUTPUT_FILEPATH',
    dest='output',
    required=True
    )

args = parser.parse_args()

def find_date_in_string(input_string, date_pattern):
    """Searches for a date within a given string."""
    date = "(No date found)"
    match = re.search(date_pattern, input_string)
    if match:
        date_matched = match.group(1)
        if len(date_matched) == 8:
            date = datetime.strptime(date_matched, "%Y%m%d").strftime("%d-%m-%Y")
        elif len(date_matched) > 8:
            date = date_matched
    return date

def get_sample_ids(samplesheet_csv):
    """Get sample id from csv."""
    df = pd.read_csv(samplesheet_csv)
    sample_ids = df['sample'].tolist()
    return sample_ids

def generate_master_html(template_html_fpath, sample_ids, seqrun_date, timestamp):
    """Read the template from an HTML file."""
    with open(template_html_fpath, "r") as file:
        master_template = file.read()
    template = Template(master_template)
    rendered_html = template.render(sample_ids=sample_ids, seqrun_date=seqrun_date, timestamp=timestamp)
    return rendered_html

def main():
    sample_ids = get_sample_ids(args.csv)
    seqrun_date = find_date_in_string(args.csv, r'/(\d{8})_')
    rendered_html = generate_master_html(args.html, sample_ids, seqrun_date, args.timestamp)
    with open(args.output, "w") as fout:
        fout.write(rendered_html)

if __name__ == "__main__":
    main()
