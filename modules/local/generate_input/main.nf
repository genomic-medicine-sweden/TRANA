

// Merge nanopore barcode fastq.gz files when you have have sample sheet for the barcode folders
process GENERATE_INPUT {
    debug true //print to stdout. debugging

    //               Software MUST be pinned to channel (i.e. "bioconda"), version (i.e. "1.10").
    //               For Conda, the build (i.e. "h9402c20_2") must be EXCLUDED to support installation on different operating systems.
    conda "conda-forge::python=3.9"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.9':
        'quay.io/biocontainers/python:3.9' }"

    input:
    path(merged_files)

    output:
    // publishDir 'fastq_pass_merged', mode: 'move'
    path '*amplesheet_merged.csv' , emit : sample_sheet_merged

    script:
    """
    generate_input.sh $merged_files
    """

    stub:
    """
    touch samplesheet_merged.csv
    """
}

