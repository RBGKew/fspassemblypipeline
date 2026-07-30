process GETSEQKITK {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/coreutils:9.5' :
        'quay.io/biocontainers/coreutils:9.5' }"

    input:
    tuple val(meta), path(stats)

    output:
    tuple val(meta), path("*.txt"), emit: seqkitkmer_txt
    tuple val("${task.process}"), val('bash'), eval("bash --version | head -1 | sed 's/.*version //; s/ .*//'"), topic: versions, emit: versions_getseqkitk

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    awk '{print \$10}' ${stats} | \\
        tail -n1 | \\
        awk '{result = int(\$1 * 2 / 3); print (result % 2 == 0) ? result - 1 : result}' \\
        > ${prefix}_seqkit_kmer.txt
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "21" > ${prefix}_seqkit_kmer.txt
    touch ${prefix}_seqkit_kmer.txt
    """
}
