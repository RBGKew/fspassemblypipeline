process EXTRACTASSEMBLYMETRICS {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/81/81fd170d9d2bd4aa6561ccd1638a79e27895e4884753200237b6695e626794d2/data':
        'community.wave.seqera.io/library/matplotlib_numpy_pandas_python:a03b8e12185b7953' }"

    input:
    tuple val(meta), path(busco_general), path(busco_specific),
          path(quast_results), path(merqury_qv), path(merqury_completeness),
          path(coverage_summary)

    output:
    path("${meta.id}_assembly_metrics.tsv"), emit: metrics_row
    tuple val("${task.process}"), val('python'), eval("python --version"), topic: versions, emit: versions_extractassemblymetrics

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    extract_metrics.py \\
        --sample-id ${meta.id} \\
        --reads-type ${meta.reads_type} \\
        --kmer-strategy ${meta.kmer_strategy} \\
        --assembler ${meta.assembler} \\
        --polisher pypolca \\
        --busco-general ${busco_general} \\
        --busco-specific ${busco_specific} \\
        --quast-tsv ${quast_results}/transposed_report.tsv \\
        --merqury-qv ${merqury_qv} \\
        --merqury-completeness ${merqury_completeness} \\
        --coverage-summary ${coverage_summary} \\
        --output ${meta.id}_assembly_metrics.tsv
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    echo $args

    touch ${prefix}_assembly_metrics.tsv
    """
}
