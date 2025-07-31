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


process ASSIGNMENT_HEATMAP {
    debug true
    tag "$meta.id"
    label 'process_single'

    //               Software MUST be pinned to channel (i.e. "bioconda"), version (i.e. "1.10").
    //               For Conda, the build (i.e. "h9402c20_2") must be EXCLUDED to support installation on different operating systems.
    conda 'modules/local/assignment_heatmap/env.yaml'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/15/152f64a06bd7bf260139b6eae0bf7c6c7bc7f3e13011d9b70a32cc03e0986250/data':
        'community.wave.seqera.io/library/bioconductor-phyloseq_r-base_r-stringr_r-tidyr:2c6c1953585e3a79' }"

    input:
    //  Where applicable all sample-specific information e.g. "id", "single_end", "read_group"
    //               MUST be provided as an input via a Groovy Map called "meta".
    //               This information may not be required in some instances e.g. indexing reference genome files:
    //               https://github.com/nf-core/modules/blob/master/modules/nf-core/bwa/index/main.nf
    //  Where applicable please provide/convert compressed files as input/output
    //               e.g. "*.fastq.gz" and NOT "*.fastq", "*.bam" and NOT "*.sam" etc.
    tuple val(meta), path(assignment_report)

    output:
    path "*assignment_heatmap.png"    , emit: assignment_heatmap
    path "versions.yml"             , emit: versions
    path "*assignment_heatmap_log.txt"             , emit: assignment_heatmap_log

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    assignment_heatmap.R  "$assignment_report" "${prefix}_assignment_heatmap.png" > ${prefix}_assignment_heatmap_log.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(echo \$(R --version 2>&1) | sed 's/^.*R version //; s/ .*\$//')
        bioconductor-phyloseq: \$(Rscript -e "library(phyloseq); cat(as.character(packageVersion('phyloseq')))")
        r-tidyr: \$(Rscript -e "library(tidyr); cat(as.character(packageVersion('tidyr')))")
        r-stringr: \$(Rscript -e "library(stringr); cat(as.character(packageVersion('stringr')))")
        r-stringr: \$(Rscript -e "library(ggplot2); cat(as.character(packageVersion('ggplot2')))")
        r-stringr: \$(Rscript -e "library(data.table); cat(as.character(packageVersion('data.table')))")

    END_VERSIONS
    """
}

