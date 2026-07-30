process COMPARISON {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "community.wave.seqera.io/library/r-base_r-data.table_r-dplr_r-stringr_r-tidyverse:27b5ede88e4d8afb"

    input:
    tuple val(meta), path(classifications), path(fcs_report_reformatted)

    output:
    tuple val(meta), path("*_tiara_vs_fcs_compare.tsv"), emit: comparison_table
    tuple val(meta), path("*_blobtools_taxonomy.tsv")  , emit: blobtools_taxonomy
    tuple val("${task.process}"), val('comparison'), eval("R --version | head -n1 | sed 's/R version //; s/ .*//'"), topic: versions, emit: versions_comparison

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    comparison.R \\
        ${fcs_report_reformatted} \\
        ${classifications} \\
        ${prefix}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_tiara_vs_fcs_compare.tsv
    touch ${prefix}_blobtools_taxonomy.tsv
    """
}
