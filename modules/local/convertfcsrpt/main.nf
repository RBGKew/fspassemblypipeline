process CONVERTFCSRPT {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/coreutils:9.5' :
        'quay.io/biocontainers/coreutils:9.5' }"

    input:
    tuple val(meta), path(rpt)

    output:
    tuple val(meta), path("*.tsv"), emit: fcs_report_reformatted
    tuple val("${task.process}"), val('convertfcsrpt'), eval("cut --version | head -n1 | sed 's/cut (GNU coreutils) //'"), topic: versions, emit: versions_convertfcsrpt

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    cut --complement -f 5,12,18,24,30 ${rpt} \\
        | tail -n +2 \\
        | sed '1s/^#//' \\
        > ${prefix}.tsv
    """

}
