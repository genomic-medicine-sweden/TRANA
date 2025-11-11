

// Merge nanopore barcode fastq.gz files when you have have sample sheet for the barcode folders
process GENERATE_INPUT {
    debug false //print to stdout. debugging

    //               Software MUST be pinned to channel (i.e. "bioconda"), version (i.e. "1.10").
    //               For Conda, the build (i.e. "h9402c20_2") must be EXCLUDED to support installation on different operating systems.
    conda 'modules/local/generate_input/env.yaml'
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.9':
        'quay.io/biocontainers/python:3.9' }"

    input:
    path(merged_files)

    output:
    // publishDir 'fastq_pass_merged', mode: 'move'
    path '*amplesheet_merged.csv'                               ,emit : sample_sheet_merged
    path "versions.yml"                                         , emit: versions

    script:
    """
    {
    generate_input.sh $merged_files
    } > generate_input_log.log 2>&1

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        generate_input: \$(echo \$(generate_input.sh version 2>&1))
    END_VERSIONS
    """

    stub:
    """
    touch samplesheet_merged.csv
    touch versions.yml
    touch generate_input_log.log
    """
}



