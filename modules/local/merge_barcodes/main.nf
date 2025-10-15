// Merge nanopore barcode fastq.gz files
process MERGE_BARCODES {
    debug false //print to stdout

    //               Software MUST be pinned to channel (i.e. "bioconda"), version (i.e. "1.10").
    //               For Conda, the build (i.e. "h9402c20_2") must be EXCLUDED to support installation on different operating systems.
    conda 'modules/local/merge_barcodes/env.yaml'
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
    path "versions.yml"                                         , emit: versions
    path "merge_barcodes_log.log"                               , emit: merge_barcodes_log


    script:
    """
    {
        merge_barcodes.sh $fastq_pass fastq_pass_merged
    } > merge_barcodes_log.log 2>&1

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        merge_barcodes.sh: \$(echo \$(merge_barcodes.sh version 2>&1))
    END_VERSIONS
    """

    stub:
    """
    mkdir fastq_pass_merged
    touch fastq_pass_merged/1.fastq.gz
    touch fastq_pass_merged/merge_barcodes_log.log
    """
}

