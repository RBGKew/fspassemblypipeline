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
    def fastq_file_list = fastq_files instanceof List ? fastq_files : [fastq_files]
    def fastq_args = fastq_file_list.collect { "\"${it}\"" }.join(' ')
    """
    for fq in ${fastq_args}; do
        zcat "\$fq" | perl ${projectDir}/bin/fq_n50.pl > "${meta.id}_\$(basename "\$fq").stats"
    done
    """

    stub:
    def fastq_file_list = fastq_files instanceof List ? fastq_files : [fastq_files]
    def stat_stubs = fastq_file_list.collect { fq ->
        def base = fq.getName()
        """
        cat <<'EOF' > "${meta.id}_${base}.stats"
Total: 0
Average: 0
EOF
        """.stripIndent().trim()
    }.join('\n')
    """
    ${stat_stubs}
    """
}
