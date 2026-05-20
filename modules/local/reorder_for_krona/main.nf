process REORDER_FOR_KRONA {
    tag "$meta.id"
    
    input:
    tuple val(meta), path(abundance_file)
    
    output:
    tuple val(meta), path("${meta.id}-reordered.tsv"), emit: reordered
    
    script:
    """
    awk 'BEGIN{FS=OFS="\\t"} {print \$2,\$9,\$8,\$7,\$6,\$5,\$4,\$3}' \
      ${abundance_file} \
      > ${meta.id}-reordered.tsv
    """
}