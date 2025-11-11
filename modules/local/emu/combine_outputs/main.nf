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

process EMU_COMBINE_OUTPUTS {
    debug false
//    tag "$meta.id"
    label 'process_single'

    //               Software MUST be pinned to channel (i.e. "bioconda"), version (i.e. "1.10").
    //               For Conda, the build (i.e. "h9402c20_2") must be EXCLUDED to support installation on different operating systems.
    conda 'modules/local/emu/combine_outputs/env.yaml'
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
    path(collected_files)

    output:
    path("collected_reports_dir/emu-combined-*.tsv"), emit: combined_report
    path("collected_reports_counts_dir/emu-combined-*counts.tsv"), emit: combined_counts_report

    path "versions.yml"                             , emit: versions
    path "emu_combine_outputs.log"                  , emit: log


    when:
    task.ext.when == null || task.ext.when

    script:
    def file_list = collected_files.join(' ')

    def args = task.ext.args ?: ''
    // Uncomment and remove _ if prefix will be used.
    //def _prefix = task.ext.prefix ?: "${meta.id}"

    """
    {
    mkdir -p collected_reports_dir
    cp ${file_list} collected_reports_dir/
    mkdir -p collected_reports_counts_dir
    cp -v ${file_list} collected_reports_counts_dir/
    emu \\
        combine-outputs \\
        './collected_reports_counts_dir/' \\
        --counts \\
        $args
    emu \\
        combine-outputs \\
        './collected_reports_dir/' \\
        $args

    } > emu_combine_outputs.log 2>1
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        emu: \$(echo \$(emu --version 2>&1) | sed 's/^.*emu //; s/Using.*\$//' )
    END_VERSIONS
    """
}

