process BLOBTOOLKIT_CREATEBLOBDIR {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"

    container "docker.io/genomehubs/blobtoolkit:4.4.6"

    input:
    tuple val(meta), path(fasta), path(bam), path(busco), path(yaml), path(index), path(taxonomy)

    output:
    tuple val(meta), path("${meta.id}"), emit: blobdir
    tuple val("${task.process}"), val('blobtoolkit'), eval("btk --version | cut -d' ' -f2 | sed 's/v//'"), topic: versions, emit: versions_blobtoolkit

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def busco_arg = busco ? "--busco ${busco}" : ""
    // we have btk do seperate create and add because this was how we got it to work. in theory you can add all inputs in 1 create block, we had issues doing this.
    """
    blobtools create \\
        --fasta ${fasta} \\
        --meta ${yaml} \\
        --threads ${task.cpus} \\
        ${busco_arg} \\
        ${args} \\
        ${prefix}

    blobtools add \\
        --cov ${bam} \\
        --threads ${task.cpus} \\
        ${prefix}

    blobtools add \\
       --text ${taxonomy} \\
       --text-delimiter '\t' \\
       --text-cols 'seq_id=identifiers,taxonomy=taxonomy' \\
       --text-header \\
        --key plot.cat=taxonomy \\
        ${prefix}
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p ${prefix}
    touch ${prefix}/meta.json
    touch ${prefix}/identifiers.json
    """
}
