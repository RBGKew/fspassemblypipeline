process FQSTAT_SUMMARY {
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/coreutils:9.5' :
        'quay.io/biocontainers/coreutils:9.5' }"

    input:
    path(stats_files)

    output:
    path('z_states_for_spreadsheet/*txt'), emit: summary

    script:
    """
    mkdir -p z_states_for_spreadsheet
    : > z_states_for_spreadsheet/total_bp_merged.txt
    : > z_states_for_spreadsheet/Len_avg_merged.txt
    : > z_states_for_spreadsheet/total_bp_trimmed.txt

    for stat in ${stats_files}; do
        stat_base=\$(basename "\$stat")

        if [[ "\$stat_base" == *merge.fq.gz.stats || "\$stat_base" == *merged.fastq.gz.stats ]]; then
            grep 'Total' "\$stat" >> z_states_for_spreadsheet/total_bp_merged.txt || true
            grep 'Average' "\$stat" >> z_states_for_spreadsheet/Len_avg_merged.txt || true
        fi

        if [[ "\$stat_base" == *trimmed*.fq.gz.stats || "\$stat_base" == *fastp.fastq.gz.stats ]]; then
            grep 'Total' "\$stat" >> z_states_for_spreadsheet/total_bp_trimmed.txt || true
        fi
    done
    """
}
