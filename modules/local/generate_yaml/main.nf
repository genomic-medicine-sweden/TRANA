process GENERATE_YAML {
    //               Software MUST be pinned to channel (i.e. "bioconda"), version (i.e. "1.10").
    //               For Conda, the build (i.e. "pyhdfd78af_1") must be EXCLUDED to support installation on different operating systems.
    conda 'modules/local/generate_yaml/env.yaml'
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://clinicalgenomicslund/eyrie-popup:0.3.0':
        'docker://clinicalgenomicslund/eyrie-popup:0.3.0' }"

    input:
    tuple val(meta), val(input_dir)

    output:
    path output         , emit: yaml
    path "versions.yml" , emit: versions

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    output = "${prefix}.yaml"
    """
    popup generate-config \\
        --analysis-output-dirpath ${input_dir} \\
        --sample-id ${prefix} \\
        --output ${output}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        eyrie-popup: \$(echo \$(popup --version 2>&1) | sed 's/^.*eyrie-popup, version //; s/ .*\$//' )
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    output = "${prefix}.yaml"
    """
    touch ${output}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        eyrie-popup: \$(echo \$(popup --version 2>&1) | sed 's/^.*eyrie-popup, version //; s/ .*\$//' )
    END_VERSIONS
    """
}
