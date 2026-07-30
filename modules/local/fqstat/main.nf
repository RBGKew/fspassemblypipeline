process FQSTAT {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/18/18a14d4248d24b68ddd07a6ba1988b6a9cc99a134ecce2a4e00521d44ebae200/data' :
        'community.wave.seqera.io/library/coreutils_gzip_perl:6929d8dfa75ef275' }"

    input:
    tuple val(meta), path(fastq_files)

    output:
    tuple val(meta), path("${meta.id}*.stats"), emit: stats
    tuple val("${task.process}"), val('perl'), eval("perl -v | grep -oP '\(v\K[0-9.]+'"), topic: versions, emit: versions_fqstat_summary

    script:
    def fastq_file_list = fastq_files instanceof List ? fastq_files : [fastq_files]
    def fastq_args = fastq_file_list.collect { "\"${it}\"" }.join(' ')
    """
    for fq in ${fastq_args}; do
        zcat "\$fq" | perl ${projectDir}/bin/fq_n50.pl > "${meta.id}_\$(basename "\$fq").stats"
    done
    """

    stub:
    def stat_stubs = [fastq_files].flatten().collect { fq ->
        "touch \"${meta.id}_${fq.getName()}.stats\""
    }.join('\n')
    """
    ${stat_stubs}
    """
}
