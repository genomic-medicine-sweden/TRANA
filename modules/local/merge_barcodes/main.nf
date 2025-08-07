// Merge nanopore barcode fastq.gz files
process MERGE_BARCODES {
    debug false //print to stdout

    //               Software MUST be pinned to channel (i.e. "bioconda"), version (i.e. "1.10").
    //               For Conda, the build (i.e. "h9402c20_2") must be EXCLUDED to support installation on different operating systems.
    conda "conda-forge::python=3.9"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.9':
        'quay.io/biocontainers/python:3.9' }"

    input:
    path fastq_pass

    output:
    // publishDir 'fastq_pass_merged', mode: 'move'
    // path('*fastq.gz') , emit : fastq_files_merged
    path('fastq_pass_merged/*fastq.gz') , emit : fastq_files_merged
    path('fastq_pass_merged')           , emit : fastq_dir_merged

    script:
    """
    merge_barcodes.sh $fastq_pass fastq_pass_merged
    """

    stub:
    """
    mkdir fastq_pass_merged
    touch fastq_pass_merged/1.fastq.gz
    """
}

