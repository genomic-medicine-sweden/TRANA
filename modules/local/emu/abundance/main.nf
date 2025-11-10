//  A module file SHOULD only define input and output files as command-line parameters.
//               All other parameters MUST be provided using the "task.ext" directive, see here:
//               https://www.nextflow.io/docs/latest/process.html#ext
//               where "task.ext" is a string.
//               Any parameters that need to be evaluated in the context of a particular sample
//               e.g. single-end/paired-end data MUST also be defined and evaluated appropriately.
//  Software that can be piped together SHOULD be added to separate module files
//               unless there is a run-time, storage advantage in implementing in this way
//               e.g. it's ok to have a single module for bwa to output BAM instead of SAM:
//                 bwa mem | samtools view -B -T ref.fasta
//  Optional inputs are not currently supported by Nextflow. However, using an empty
//               list (`[]`) instead of a file can be used to work around this issue.

process EMU_ABUNDANCE {
    tag "$meta.id"
    label 'process_high'

    //               Software MUST be pinned to channel (i.e. "bioconda"), version (i.e. "1.10").
    //               For Conda, the build (i.e. "h9402c20_2") must be EXCLUDED to support installation on different operating systems.
    conda 'modules/local/emu/abundance/env.yaml'
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/89/89f51e5b1961b27efc33dbafe25ae7f85c1ccfc2e2df5341849237a6db2023a1/data':
        'community.wave.seqera.io/library/emu:3.5.4--c9ac9572d0d77ae9' }"

    input:
    //  Where applicable all sample-specific information e.g. "id", "single_end", "read_group"
    //               MUST be provided as an input via a Groovy Map called "meta".
    //               This information may not be required in some instances e.g. indexing reference genome files:
    //               https://github.com/nf-core/modules/blob/master/modules/nf-core/bwa/index/main.nf
    //  Where applicable please provide/convert compressed files as input/output
    //               e.g. "*.fastq.gz" and NOT "*.fastq", "*.bam" and NOT "*.sam" etc.
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*abundance.tsv")                     , emit: report
    tuple val(meta), path("*read-assignment-distributions.tsv") , emit: assignment_report, optional:true
    tuple val(meta), path("*.sam")                              , emit: samfile, optional:true
    tuple val(meta), path("*.fastq_unclassified_mapped.fasta")  , emit: unclassified_mapped_fa , optional:true
    tuple val(meta), path("*.fastq_unmapped.fasta")             , emit: unclassified_unmapped_fa , optional:true
    path "versions.yml"                                         , emit: versions
    path "*_emu_log.log"                                        , emit: emu_log


    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    {
    echo "$reads"
    emu \\
        abundance \\
        $args \\
        --threads $task.cpus \\
        $reads
    } > "${prefix}_emu_log.log" 2>&1

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        emu: \$(echo \$(emu --version 2>&1) | sed 's/^.*emu //; s/Using.*\$//' )
    END_VERSIONS
    """

    stub:
    def prefix  = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.abundance.tsv
    touch ${prefix}.read-assignment-distributions.tsv
    touch ${prefix}.sam
    touch ${prefix}.fastq_unclassified_mapped.fasta
    touch ${prefix}.fastq_unmapped.fasta
    touch ${prefix}_emu_log.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        emu: \$(echo \$(emu --version 2>&1) | sed 's/^.*emu //; s/Using.*\$//' )
    END_VERSIONS
    """
}

