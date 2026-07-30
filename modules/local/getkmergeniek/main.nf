process GETKMERGENIEK {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/coreutils:9.5' :
        'quay.io/biocontainers/coreutils:9.5' }"

    input:
    tuple val(meta), path(html)

    output:
    tuple val(meta), path("*.txt"), emit: kmer_txt
    tuple val("${task.process}"), val('bash'), eval("bash --version | head -1 | sed 's/.*version //; s/ .*//'"), topic: versions, emit: versions_getkmergeniek

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    extract_best_k.awk < ${html} > ${prefix}_best_kmer.txt
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "21" > ${prefix}_best_kmer.txt
    touch ${prefix}_best_kmer.txt
    """
}
