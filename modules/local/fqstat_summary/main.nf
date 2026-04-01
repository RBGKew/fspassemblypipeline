process FQSTAT_SUMMARY {
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/coreutils:9.5' :
        'quay.io/biocontainers/coreutils:9.5' }"

    input:
    path(stats_files)

    output:
    path('fq_stats.summary.*.txt'), emit: summary

    script:
    def stats_file_list = stats_files instanceof List ? stats_files : [stats_files]
    def stats_args = stats_file_list.collect { "\"${it}\"" }.join(' ')
    def merged_summary_file = 'fq_stats.summary.merged.txt'
    def trimmed_summary_file = 'fq_stats.summary.trimmed.txt'
    """
    printf "stats_files\ttotal_bp\tLen_avg\n" > ${merged_summary_file}
    printf "stats_files\ttotal_bp\tLen_avg\n" > ${trimmed_summary_file}

    for stat in ${stats_args}; do
        stat_base=\$(basename "\$stat")
        total_bp=\$(awk '\$1=="Total:" {print \$2; exit}' "\$stat")
        Len_avg=\$(awk '\$1=="Average:" {print \$2; exit}' "\$stat")

        if [[ "\$stat_base" == *merge.fq.gz.stats || "\$stat_base" == *merged.fastq.gz.stats ]]; then
            printf "%s\t%s\t%s\n" "\$stat_base" "\$total_bp" "\$Len_avg" >> ${merged_summary_file}
        fi

        if [[ "\$stat_base" == *trimmed*.fq.gz.stats || "\$stat_base" == *fastp.fastq.gz.stats ]]; then
            printf "%s\t%s\t%s\n" "\$stat_base" "\$total_bp" "\$Len_avg" >> ${trimmed_summary_file}
        fi
    done
    """
}
