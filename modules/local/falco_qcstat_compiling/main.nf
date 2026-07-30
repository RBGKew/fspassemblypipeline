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
    tuple val(qc_subdir), path('QC_*_result.txt'), emit: stats
    tuple val("${task.process}"), val('bash'), eval("bash --version | head -1 | sed 's/.*version //; s/ .*//'"), topic: versions, emit: versions_falco_qcstat_compiling

    script:
    def falco_txt_list = falco_txt_files instanceof List ? falco_txt_files : [falco_txt_files]
    def falco_txt_args = falco_txt_list
        .collect { "\"${it}\"" }
        .join(' ')
    def output_name_map = [
        raw: 'QC_raw_result.txt',
        trimmed: 'QC_trimmed_result.txt',
        merge: 'QC_merge_result.txt'
    ]
    def output_name = output_name_map[qc_subdir as String] ?: "QC_${qc_subdir}_result.txt"
    """
    bash ${projectDir}/bin/fastQC_result_compiling.sh ${falco_txt_args}
    mv fastQC_result.txt ${output_name}
    """
}
