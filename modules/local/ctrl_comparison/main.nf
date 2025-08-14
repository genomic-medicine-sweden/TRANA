process CTRL_COMPARISON {
    debug false
    //tag "$meta.id"
    label 'process_single'

    //               Software MUST be pinned to channel (i.e. "bioconda"), version (i.e. "1.10").
    //               For Conda, the build (i.e. "h9402c20_2") must be EXCLUDED to support installation on different operating systems.
    conda 'modules/local/assignment_heatmap/env.yaml'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/70/701063393e92e8ae9c7632857dd6f5c42c04f6d7f26f7ea299cb89c98eae9226/data':
        'community.wave.seqera.io/library/r-base_r-dplyr_r-ggplot2_r-glue_pruned:c573c918fb79c025'}"

    input:
    //  Where applicable all sample-specific information e.g. "id", "single_end", "read_group"
    //               MUST be provided as an input via a Groovy Map called "meta".
    //               This information may not be required in some instances e.g. indexing reference genome files:
    //               https://github.com/nf-core/modules/blob/master/modules/nf-core/bwa/index/main.nf
    //  Where applicable please provide/convert compressed files as input/output
    //               e.g. "*.fastq.gz" and NOT "*.fastq", "*.bam" and NOT "*.sam" etc.
    path(combined_report)
    path(combined_counts_report)


    output:
    path "plots_vs_selected_ctrl/*.png"    , emit: ctrl_comparison_png
    path "versions.yml"                    , emit: versions
    path "*ctrl_comparison_log.log"        , emit: ctrl_comparison_log

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def _prefix = task.ext.prefix ?: "${meta.id}"
    """
    {
    export XDG_CACHE_HOME=".cache"
    mkdir -p \${XDG_CACHE_HOME}/fontconfig

    # ctrl_comparison_cli.R ctrl_comparison_cli_data.tsv my_neg.fastq my_pos.fastq
    ctrl_comparison_cli.R \\
    $combined_report \\
    $args

    ctrl_comparison_cli.R \\
    $combined_counts_report \\
    $args \\
    --counts_file

    } > ctrl_comparison_log.log 2>&1

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ctrl_comparison_cli.R: \$(echo \$(ctrl_comparison_cli.R --version 2>&1)  | sed 's/ctrl_comparison_cli.R version //')
        r-ggplot2: \$(Rscript -e "library(ggplot2); cat(as.character(packageVersion('ggplot2')))")
        r-dplyr: \$(Rscript -e "library(dplyr); cat(as.character(packageVersion('dplyr')))")
        r-glue: \$(Rscript -e "library(glue); cat(as.character(packageVersion('glue')))")
        r-readr: \$(Rscript -e "library(readr); cat(as.character(packageVersion('readr')))")
        r-tidyr: \$(Rscript -e "library(tidyr); cat(as.character(packageVersion('tidyr')))")
        r-viridis: \$(Rscript -e "library(viridis); cat(as.character(packageVersion('viridis')))")
        r-viridisLite: \$(Rscript -e "library(viridisLite); cat(as.character(packageVersion('viridisLite')))")
    END_VERSIONS
    """
}

