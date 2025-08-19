process TRANSLATE_TAXIDS {
    debug false
    tag "$meta.id"
    label 'process_single'

    //               Software MUST be pinned to channel (i.e. "bioconda"), version (i.e. "1.10").
    //               For Conda, the build (i.e. "h9402c20_2") must be EXCLUDED to support installation on different operating systems.
    conda 'modules/local/translate_taxids/env.yaml'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/63/6325e0e1cde97907c30be26d33a170477cca271da36267bf0679965d1b08f0d2/data':
        'community.wave.seqera.io/library/pandas_python:4330fd07d14e9bfb' }"

    input:
    //  Where applicable all sample-specific information e.g. "id", "single_end", "read_group"
    //               MUST be provided as an input via a Groovy Map called "meta".
    //               This information may not be required in some instances e.g. indexing reference genome files:
    //               https://github.com/nf-core/modules/blob/master/modules/nf-core/bwa/index/main.nf
    //  Where applicable please provide/convert compressed files as input/output
    //               e.g. "*.fastq.gz" and NOT "*.fastq", "*.bam" and NOT "*.sam" etc.
    tuple val(meta), path(assignment_report)


    output:
    tuple val(meta), path("*read-assignment-distributions_translated.tsv") , emit: assignment_translated_report
    path "versions.yml"                  , emit: versions
    path "*translate_taxids_log.log"   , emit: translate_taxids_log

    when:
    task.ext.when == null || task.ext.when

    script:
    def _args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    {
    translate_taxids.py $assignment_report  "${params.db}/taxonomy.tsv" ${prefix}_read-assignment-distributions_translated.tsv
    } > translate_taxids_log.log 2>&1
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        translate_taxids: \$(translate_taxids.py --version | sed 's/translate_taxids.py version //')
    END_VERSIONS
    """
}



