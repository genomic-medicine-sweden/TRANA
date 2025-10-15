// Merge nanopore barcode fastq.gz files when you have have sample sheet for the barcode folders
process MERGE_BARCODES_SAMPLESHEET {
    debug false //print to stdout. debugging

    //               Software MUST be pinned to channel (i.e. "bioconda"), version (i.e. "1.10").
    //               For Conda, the build (i.e. "h9402c20_2") must be EXCLUDED to support installation on different operating systems.
    conda 'modules/local/merge_barcodes_samplesheet/env.yaml'
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/nf-core:3.0.2--pyhdfd78af_1':
        'quay.io/biocontainers/nf-core:3.0.2--pyhdfd78af_1' }"

    input:
    path barcodes_samplesheet
    path fastq_pass

    output:
    // publishDir 'fastq_pass_merged'   , mode: 'move'
    path('fastq_pass_merged/*fastq.gz') , emit : fastq_files_merged
    path('fastq_pass_merged')           , emit : fastq_dir_merged
    path "versions.yml"                 , emit: versions

    script:
    """
    merge_barcodes_samplesheet.py ${barcodes_samplesheet} fastq_pass_merged ${fastq_pass}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        merge_barcodes_samplesheet.py: \$(echo \$(merge_barcodes_samplesheet.py --version 2>&1) | sed 's/^.*merge_barcodes_samplesheet.py, version //; s/ .*\$//' )
    END_VERSIONS
    """

    stub:
    """
    mkdir fastq_pass_merged
    touch fastq_pass_merged/1.fastq.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        merge_barcodes_samplesheet.py: \$(echo \$(merge_barcodes_samplesheet.py --version 2>&1) | sed 's/^.*merge_barcodes_samplesheet.py, version //; s/ .*\$//' )
    END_VERSIONS
    """
}

