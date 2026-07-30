process CREATE_PROJECT_YAML {
    tag "$meta.id"
    label 'process_single'

    conda "conda-forge::python=3.11"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/coreutils:9.5' :
        'quay.io/biocontainers/coreutils:9.5' }"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("${meta.id}.yaml"), emit: yaml
    tuple val("${task.process}"), val('bash'), eval("bash --version | head -1 | sed 's/.*version //; s/ .*//'"), topic: versions, emit: versions_create_project_yaml

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def species = meta.species ?: "Unknown species"
    def taxon_id = meta.taxid ?: 1
    def accession = meta.accession ?: meta.id
    def prefix = meta.prefix ?: meta.id
    def record_type = meta.record_type ?: "contig"
    """
    cat > ${meta.id}.yaml <<EOF
assembly:
  accession: ${accession}
  prefix: ${prefix}
  record_type: ${record_type}
taxon:
  name: ${species}
  taxid: ${taxon_id}
EOF
    """

    stub:
    def species = meta.species ?: "Unknown species"
    def taxon_id = meta.taxid ?: 1
    def accession = meta.accession ?: meta.id
    def prefix = meta.prefix ?: meta.id
    def record_type = meta.record_type ?: "contig"
    """
    cat > ${meta.id}.yaml <<EOF
assembly:
  accession: ${accession}
  prefix: ${prefix}
  record_type: ${record_type}
taxon:
  name: ${species}
  taxid: ${taxon_id}
EOF
    """
}
