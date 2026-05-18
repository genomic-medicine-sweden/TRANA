process REORDER_FOR_KRONA {
    tag "$meta.id"
    
    input:
    tuple val(meta), path(abundance_file)
    
    output:
    tuple val(meta), path("${meta.id}-reordered.tsv"), emit: reordered
    
    script:
    """
    cat ${abundance_file} \\
      | cut -f 2,8,7,6,5,4,3 \\
      > ${meta.id}-reordered.tsv
    """
}