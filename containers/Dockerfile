FROM python:3.11-slim

COPY bin/check_samplesheet.py bin/generate_master_html.py bin/merge_barcodes_samplesheet.py bin/requirements*.txt .

# install dependancies
RUN pip install --no-cache-dir -r requirements.txt
