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
process COMBINE_REPORTS {
    label 'process_single'

    input:
    //  Where applicable all sample-specific information e.g. "id", "single_end", "read_group"
    //               MUST be provided as an input via a Groovy Map called "meta".
    //               This information may not be required in some instances e.g. indexing reference genome files:
    //               https://github.com/nf-core/modules/blob/master/modules/nf-core/bwa/index/main.nf
    //  Where applicable please provide/convert compressed files as input/output
    //               e.g. "*.fastq.gz" and NOT "*.fastq", "*.bam" and NOT "*.sam" etc.
    //tuple val(meta), path(report)
    path report
    // collect all reports

    output:
    path "combined-rel-abundance.tsv", emit: combinedreport


    script:
    """
    {
        bash combine_tsv.sh "${report}" > "combined-rel-abundance.tsv"
    } > combine_reports_log.log 2>1
    """
}


process PHYLOSEQ_OBJECT {
    label 'process_single'

    //               Software MUST be pinned to channel (i.e. "bioconda"), version (i.e. "1.10").
    //               For Conda, the build (i.e. "h9402c20_2") must be EXCLUDED to support installation on different operating systems.
    conda 'modules/local/phyloseq/env.yaml'

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
    path combined_report
    path taxonomy_file

    output:
    path "phyloseq_object.RData"          , emit: phyloseq_output_RData
    path "versions.yml"                   , emit: versions
    path "phyloseq_object_log.log"        , emit: phyloseq_object_log

    script:
    """
    {
        phyloseq_object.R  $combined_report $taxonomy_file
    } > phyloseq_object_log.log 2>&1
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        phyloseq_object.R: \$(echo \$(phyloseq_object.R --version 2>&1) | grep -i 'version' | sed 's/phyloseq_object.R version//')
        combine_reports.sh: \$(echo 1.0)
        r-base: \$(echo \$(R --version 2>&1) | sed 's/^.*R version //; s/ .*\$//')
        bioconductor-phyloseq: \$(Rscript -e "library(phyloseq); cat(as.character(packageVersion('phyloseq')))")
        r-tidyr: \$(Rscript -e "library(tidyr); cat(as.character(packageVersion('tidyr')))")
        r-stringr: \$(Rscript -e "library(stringr); cat(as.character(packageVersion('stringr')))")
    END_VERSIONS
    """
}

