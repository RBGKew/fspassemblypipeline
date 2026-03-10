process FALCO_QCSTAT_COMPILING {
    tag "$qc_subdir"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/coreutils:9.5' :
        'quay.io/biocontainers/coreutils:9.5' }"

    input:
    tuple val(qc_subdir), path(falco_txt_files)

    output:
    tuple val(qc_subdir), path('fastQC_result.txt'), emit: stats

    script:
    """
    mkdir -p falco_inputs
    cp ${falco_txt_files} falco_inputs/

    # The legacy compiler script expects files matching *gz_summary.txt.
    for summary in falco_inputs/*_summary.txt; do
        [[ -e "\$summary" ]] || continue
        if [[ "\$summary" != *gz_summary.txt ]]; then
            cp "\$summary" "\${summary%_summary.txt}.fq.gz_summary.txt"
        fi
    done

    bash ${projectDir}/bin/fastQC_result_compiling.sh
    """
}
