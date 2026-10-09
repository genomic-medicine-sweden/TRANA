process EMUSE {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
?         'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/af/afec8be6a94d99b6e61859c447321ff68a1006904fb4a1e3f2f54117d29ab970/data'
:         'community.wave.seqera.io/library/emuse:1.1.0--a094b001db16b1b3' }"

    input:
    tuple val(meta), path(input_dir)

    output:
    tuple val(meta), path("*.html"), emit: html
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    emuse \\
        ${args} \\
        --input-dir "${input_dir}" \\
        --output-file "${prefix}.html" \\
        --sample-name "${meta.id}" \\
        --neg-control "${meta.neg_control}"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        emuse: \$(emuse --version | sed 's/^.* //; s/^v//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.html

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        emuse: \$(emuse --version | sed 's/^.* //; s/^v//')
    END_VERSIONS
    """
}
