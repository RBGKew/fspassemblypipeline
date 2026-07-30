process COVERAGEVIZ {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/81/81fd170d9d2bd4aa6561ccd1638a79e27895e4884753200237b6695e626794d2/data':
        'community.wave.seqera.io/library/matplotlib_numpy_pandas_python:a03b8e12185b7953' }"

    input:
    tuple val(meta), path(samtools_coverage_output), path(samtools_flagstats_output)

    output:
    tuple val(meta), path("*_coverage_plot.png"),    emit: png
    tuple val(meta), path("*_coverage_summary.txt"), emit: summary
    tuple val(meta), path("${meta.id}-coverageviz.log"),      emit: log
    tuple val("${task.process}"), val('python'), eval("python --version"), topic: versions, emit: versions_coverageviz

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    coverage_visualisation.py \\
        ${samtools_coverage_output} \\
        ${prefix}_coverage_plot.png \\
        --sample ${meta.id} \\
        --assembler "reads type: ${meta.reads_type}, kmer strategy: ${meta.kmer_strategy}, assembler: ${meta.assembler}, polisher: pypolca" \\
        --flagstat ${samtools_flagstats_output} \\
        --summary ${prefix}_coverage_summary.txt > ${prefix}-coverageviz.log 2>&1
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_coverage_plot.png
    touch ${prefix}_coverage_summary.txt
    touch ${prefix}-coverageviz.log
    """
}
