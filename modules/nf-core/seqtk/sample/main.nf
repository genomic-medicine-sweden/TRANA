process SEQTK_SAMPLE {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/seqtk:1.4--he4a0461_1' :
        'quay.io/biocontainers/seqtk:1.4--he4a0461_1' }"

    input:
    tuple val(meta), path(reads), val(downsample_n_reads)

    output:
    tuple val(meta), path("*downsampled*.fastq.gz"), emit: reads
    path "versions.yml"                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    if (!(args ==~ /.*\ -s\ ?[0-9]+.*/)) {
        args += " -s100"
    }
    if ( !downsample_n_reads ) {
        error "SEQTK/SAMPLE must have a downsample_n_reads value included"
    }
    """
    printf "%s\\n" $reads | while read f;
    do
        seqtk \\
            sample \\
            $args \\
            \$f \\
            $downsample_n_reads \\
            | gzip --no-name > ${prefix}.fastq.gz
    done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        seqtk: \$(echo \$(seqtk 2>&1) | sed 's/^.*Version: //; s/ .*\$//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    echo "" | gzip > ${prefix}.fastq.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        seqtk: \$(echo \$(seqtk 2>&1) | sed 's/^.*Version: //; s/ .*\$//')
    END_VERSIONS
    """

}
