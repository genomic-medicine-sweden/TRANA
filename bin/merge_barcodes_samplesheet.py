#!/usr/bin/env python
"""
Merge fastq.gz files present in barcode folders.

Author(s): Frans Wallin and Markus H Johansson
Date: 20230707
"""
import logging
import shutil
from csv import DictReader
from io import TextIOWrapper
from pathlib import Path

import click
from pydantic import BaseModel, ValidationError

LOG = logging.getLogger(__name__)


class SampleObject(BaseModel):
    """Barcode."""

    barcode: str
    sample_id: str


@click.command()
@click.argument("sample-sheet", type=click.File())
@click.argument("output-dir", type=click.Path(dir_okay=True, readable=True, path_type=Path))
@click.argument("fastq-dir", type=click.Path(dir_okay=True, readable=True, path_type=Path))
@click.version_option(version="0.2.0", prog_name="merge_barcodes_samplesheet.py")
def cli(sample_sheet: TextIOWrapper, output_dir: Path, fastq_dir: Path) -> None:
    """Merge fastq files on barcodes in a sample sheet."""
    # Create output folder and make ensure its writeable
    try:
        output_dir.mkdir(parents=True, exist_ok=True)
    except PermissionError:
        raise click.UsageError(f"Could not create '{output_dir}', output directory is not writable.")

    # Read the sample sheet file and validate the contents
    click.secho(f"Using sample sheet: {sample_sheet.name}")
    csv_file = DictReader(sample_sheet, delimiter=",", fieldnames=["barcode", "sample_id"])
    samples: list[SampleObject] = []
    for row_no, row in enumerate(csv_file, start=1):
        # skip suspected header
        if row["barcode"] == "barcode":
            continue

        # validate file
        try:
            sample_obj = SampleObject(barcode=row["barcode"], sample_id=row["sample_id"])
            # assert that barcode is in fastq directory
            if not fastq_dir.joinpath(sample_obj.barcode).exists():
                raise FileNotFoundError(sample_obj.barcode)
        except ValidationError:
            click.Abort(f"Malformed barcode or sample id at line {row_no} in {sample_sheet.name}")
        except FileNotFoundError:
            LOG.error("Barcode '%s' not found in fastq directory '%s'", sample_obj.barcode, fastq_dir.absolute())
            continue
        else:
            samples.append(sample_obj)

    # merge split fastq into a single file per barcode
    for sample_obj in samples:
        merge_fastq(
            barcode_dir=fastq_dir.joinpath(sample_obj.barcode),
            output_file=output_dir.joinpath(f"{sample_obj.sample_id}.fastq.gz"),
        )
    click.secho(f"Merged {len(samples)} ", bg="green")


def merge_fastq(barcode_dir: Path, output_file: Path) -> Path:
    """Find split gzipped fastq files in path and merge into the output file."""
    LOG.info("Merging fq files into %s", output_file)
    with output_file.open("wb") as out_fq:
        # lookup split fq
        split_fq_files = list(barcode_dir.glob("*.fastq.gz"))
        if len(split_fq_files) == 0:
            LOG.warning("No fastq files found for barcode '%s'", barcode_dir.name)

        # merge files
        for fq_file in split_fq_files:
            with fq_file.open("rb") as fq:
                LOG.debug("Merging %s file into %s", fq_file, output_file)
                shutil.copyfileobj(fq, out_fq)
    return output_file


if __name__ == "__main__":
    cli()
