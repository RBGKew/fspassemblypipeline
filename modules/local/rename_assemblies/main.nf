// I need to rename the assemblies and standardise their format to gzipped.

process RENAME_ASSEMBLIES {
    tag "$meta.id"
    label 'process_low'
    
    conda "${moduleDir}/environment.yml"
//    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
//        'https://depot.galaxyproject.org/singularity/my_tool:1.0--hdfd78af_0' :
//        'quay.io/biocontainers/my_tool:1.0--hdfd78af_0' }"

    input:
    tuple val(meta), path(input_file), val(new_name)

    output:
    tuple val(meta), path("${new_name}.gz"), emit: renamed_assemblies

    script:
    def is_gzipped = input_file.toString().endsWith('.gz')
    """
    if [ "${is_gzipped}" == "true" ]; then
        # Already gzipped, just rename
        ln -s ${input_file} ${new_name}.gz
    else
        # Gzip and rename
        gzip -c ${input_file} > ${new_name}.gz
    fi
    """
}