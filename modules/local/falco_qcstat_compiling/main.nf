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
    def falco_txt_list = falco_txt_files instanceof List ? falco_txt_files : [falco_txt_files]
    def falco_txt_args = falco_txt_list
        .collect { "\"${it}\"" }
        .join(' ')
    """
    bash ${projectDir}/bin/fastQC_result_compiling.sh ${falco_txt_args}
    """
}
