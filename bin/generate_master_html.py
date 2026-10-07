#!/usr/bin/env python

"""Generate a master html template."""

import re
import argparse
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
    find_date_in_string, generate_master_html.

List of standard modules:
    re, argparse.

List of "non standard" modules:
    jinja2.

Procedure:
    1. Get sample IDs passed in from the workflow.
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
    '-s', '--sample-ids',
    help='comma-separated list of sample ids',
    metavar='SAMPLE_IDS',
    dest='sample_ids',
    required=True
    )
parser.add_argument(
    '-p', '--path',
    help='representative input path used to detect the sequencing run date',
    metavar='PATH',
    dest='path',
    default=''
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

def generate_master_html(template_html_fpath, sample_ids, seqrun_date, timestamp):
    """Read the template from an HTML file."""
    with open(template_html_fpath, "r") as file:
        master_template = file.read()
    template = Template(master_template)
    rendered_html = template.render(sample_ids=sample_ids, seqrun_date=seqrun_date, timestamp=timestamp)
    return rendered_html

def main():
    sample_ids = args.sample_ids.split(',')
    seqrun_date = find_date_in_string(args.path, r'/(\d{8})_') if args.path else "(No date found)"
    rendered_html = generate_master_html(args.html, sample_ids, seqrun_date, args.timestamp)
    with open(args.output, "w") as fout:
        fout.write(rendered_html)

if __name__ == "__main__":
    main()
