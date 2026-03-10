process FQSTAT {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/perl:5.32.1--9e3c43247be68b3b' :
        'community.wave.seqera.io/library/perl:5.32.1--a61125adac4a9f65' }"

    input:
    tuple val(meta), path(fastq_files)

    output:
    tuple val(meta), path("${meta.id}*.stats"), emit: stats

    script:
    """
    for fq in ${fastq_files}; do
        perl ${projectDir}/bin/fq_n50.pl "\$fq" > "${meta.id}_\$(basename "\$fq").stats"
    done
    """
}
