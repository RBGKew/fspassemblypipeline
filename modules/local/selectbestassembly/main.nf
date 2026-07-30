process SELECTBESTASSEMBLY {

    tag "$meta.id"
    label 'process_single'

    conda "conda-forge::bash=5.2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/coreutils:9.5' :
        'quay.io/biocontainers/coreutils:9.5' }"

    input:
    tuple val(meta), path(busco_summaries), path(quast_tsv), path(assemblies)

    output:
    tuple val(meta), path("${meta.id}_best_assembly.fa.gz"),           emit: best_assembly
    tuple val(meta), path("best_assembly.txt"),                        emit: best_assembly_label
    tuple val(meta), path("best_assembly_meta.txt"),                   emit: best_assembly_meta
    tuple val(meta), path("complete_single_copy_buscos.txt"),          emit: busco_scores
    tuple val(meta), path("auN_quast.txt"),                            emit: aun_scores, optional: true
    tuple val(meta), path("${meta.id}-select_best_assembly.log"),      emit: selectbestassembly_log
    tuple val("${task.process}"), val('bash'), eval("bash --version | head -1 | sed 's/.*version //; s/ .*//'"), topic: versions, emit: versions_selectbestassembly

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def busco_files_str = busco_summaries instanceof List ? busco_summaries.join(' ') : busco_summaries
    def assembly_str    = assemblies instanceof List      ? assemblies.join(' ')      : assemblies
    """
    select_best_assembly.sh \\
        -s ${prefix} \\
        -b "${busco_files_str}" \\
        -q ${quast_tsv} \\
        -a "${assembly_str}" \\
        -o . > ${prefix}-select_best_assembly.log 2>&1
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_best_assembly.fa.gz
    echo "R1R2_kmergenie_spades" > best_assembly.txt
    printf "reads_type=R1R2\nkmer_strategy=kmergenie\nassembler=spades\n" > best_assembly_meta.txt
    echo "R1R2_kmergenie_spades 100" > complete_single_copy_buscos.txt
    touch ${meta.id}-select_best_assembly.log
    """
}
