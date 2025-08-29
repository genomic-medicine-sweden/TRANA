process ASSIGNMENT_HEATMAP {
    debug false
    tag "$meta.id"
    label 'process_single'

    //               Software MUST be pinned to channel (i.e. "bioconda"), version (i.e. "1.10").
    //               For Conda, the build (i.e. "h9402c20_2") must be EXCLUDED to support installation on different operating systems.
    conda 'modules/local/assignment_heatmap/env.yaml'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/35/358914ad2acb00157affa657136b76a5eb3a196d9d2236ca215c70919e917d03/data':
        'community.wave.seqera.io/library/r-base_r-data.table_r-getopt_r-ggplot2:4290b80604e6066a'}"

    input:
    //  Where applicable all sample-specific information e.g. "id", "single_end", "read_group"
    //               MUST be provided as an input via a Groovy Map called "meta".
    //               This information may not be required in some instances e.g. indexing reference genome files:
    //               https://github.com/nf-core/modules/blob/master/modules/nf-core/bwa/index/main.nf
    //  Where applicable please provide/convert compressed files as input/output
    //               e.g. "*.fastq.gz" and NOT "*.fastq", "*.bam" and NOT "*.sam" etc.
    tuple val(meta), path(assignment_translated_report)

    output:
    path "*assignment_heatmap.png"       , emit: assignment_heatmap
    path "versions.yml"                  , emit: versions
    path "*assignment_heatmap_log.log"   , emit: assignment_heatmap_log

    when:
    task.ext.when == null || task.ext.when

    script:
    // uncomment and remove _ if you want to use args
    //def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    export XDG_CACHE_HOME=".cache"
    mkdir -p \${XDG_CACHE_HOME}/fontconfig
    assignment_heatmap.R  --input_file "$assignment_translated_report" --output_file "${prefix}_assignment_heatmap.png" > ${prefix}_assignment_heatmap_log.log 2>&1

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        assignment_heatmap: \$(assignment_heatmap.R --version | sed 's/assignment_heatmap.R version//')
        r-base: \$(echo \$(R --version 2>&1) | sed 's/^.*R version //; s/ .*\$//')
        r-ggplot2: \$(Rscript -e "library(ggplot2); cat(as.character(packageVersion('ggplot2')))")
        r-data.table: \$(Rscript -e "library(data.table); cat(as.character(packageVersion('data.table')))")
        r-getopt: \$(Rscript -e "library(getopt); cat(as.character(packageVersion('getopt')))")
    END_VERSIONS
    """
}

