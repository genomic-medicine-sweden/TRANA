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
    tuple val("${task.process}"), val('emuse'), eval("emuse --version"), topic: versions, emit: versions_emuse

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    // TODO nf-core: Where possible, a command MUST be provided to obtain the version number of the software e.g. 1.10
    //               If the software is unable to output a version number on the command-line then it can be manually specified
        --input-dir "$input_dir" \\
        --output-file "${prefix}.report" \\
        --sample-name "${meta.id}" \\
        --neg-control "${meta.neg_control}"
    //               using the Nextflow "task" variable e.g. "--threads $task.cpus"
    // TODO nf-core: Please replace the example samtools command below with your module's command
    // TODO nf-core: Please indent the command appropriately (4 spaces!!) to help with readability ;)
    emuse \\

    touch "${prefix}.report"
    // TODO nf-core: A stub section should mimic the execution of the original module as best as possible
    //               Have a look at the following examples:
    //               Simple example: https://github.com/nf-core/modules/blob/624977dfaf562211e68a8a868ca80acc8461f1ac/modules/nf-core/cutadapt/main.nf#L34-L46
    //               Complex example: https://github.com/nf-core/modules/blob/88d43dad73a675e66bff49ebb57fe657a5909018/modules/nf-core/bedtools/split/main.nf#L32-L43
    // TODO nf-core: If the module doesn't use arguments ($args), you SHOULD remove:
    //               - The definition of args `def args = task.ext.args ?: ''` above.
    //               - The use of the variable in the script `echo $args ` below.
    """
    echo $args

    touch ${prefix}.bam
    """
}
