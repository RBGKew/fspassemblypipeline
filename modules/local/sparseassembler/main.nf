process SPARSEASSEMBLER {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/sparseassembler:20160205--h9948957_10':
        'biocontainers/sparseassembler:20160205--h9948957_10' }"

    input:
    tuple val(meta), path(reads)
    val kmersize
    val genome_size
    val expected_coverage

    output:
    tuple val(meta), path("*.contigs.fa.gz"), emit: contigs
    tuple val(meta), path("*.scaffolds.fa.gz"), optional: true, emit: scaffolds
    tuple val(meta), path("*-sparseassembler.log"), emit: log
    tuple val("${task.process}"), val('sparseassembler'), eval("echo '20160205'"), topic: versions, emit: versions_sparseassembler

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def input_reads = ""
    if (meta.single_end) {
        input_reads = "f ${prefix}_input.fastq"
    } else {
        input_reads = "p1 ${prefix}_R1.fastq p2 ${prefix}_R2.fastq"
    }
    """
    # Decompress input reads
    decompress() { [[ "\$1" == *.gz ]] && zcat "\$1" || cat "\$1"; }

    if [ "${meta.single_end}" = "true" ]; then
    decompress "${reads}" > ${prefix}_input.fastq
    else
    decompress "${reads[0]}"  > ${prefix}_R1.fastq
    decompress "${reads[1]}"  > ${prefix}_R2.fastq
    fi

    SparseAssembler \\
        $args \\
        k $kmersize \\
        GS $genome_size \\
        $input_reads \\
        ExpCov $expected_coverage > ${prefix}-sparseassembler.log 2>&1

    # Compress outputs if they exist
    if [ -f Contigs.txt ]; then
        mv Contigs.txt ${prefix}.contigs.fa
        gzip -c ${prefix}.contigs.fa > ${prefix}.contigs.fa.gz
        rm ${prefix}.contigs.fa
    fi

    if [ -f SuperContigs.txt ]; then
        mv SuperContigs.txt ${prefix}.scaffolds.fa
        gzip -c ${prefix}.scaffolds.fa > ${prefix}.scaffolds.fa.gz
        rm ${prefix}.scaffolds.fa
    fi
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo $args

    touch ${prefix}.contigs.fa.gz
    touch ${prefix}.scaffolds.fa.gz
    touch ${prefix}-sparseassembler.log
    """
}
