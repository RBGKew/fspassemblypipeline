/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { GENOME_ASSEMBLY        } from '../subworkflows/local/genome_assembly/main'
include { PREPROCESSING          } from '../subworkflows/local/preprocessing/main'
include { CONTAMINATION_DETECTION} from '../subworkflows/local/contamination_detection/main'
include { MULTIQC                } from '../modules/nf-core/multiqc/main'
include { paramsSummaryMap       } from 'plugin/nf-schema'
include { paramsSummaryMultiqc   } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_fspassemblypipeline_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow FSPASSEMBLYPIPELINE {

    take:
    ch_samplesheet // channel: reads for preprocessing/assembly [[meta], [fastq files]]
    ch_fasta       // channel: genome assemblies for contamination detection [[meta], fasta]
    ch_bam         // channel: BAM files from workflow emits or other sources [[meta], bam]

    main:

    ch_versions = channel.empty()
    ch_multiqc_files = channel.empty()

    //
    // SUBWORKFLOW: Run PREPROCESSING on raw reads only
    //
    // Branch samplesheet by type (raw vs cleaned vs bam)
    ch_samplesheet
        .branch { meta, files ->
            raw: meta.type == 'raw'
                return [meta, files]
            cleaned: meta.type == 'cleaned'
                return [meta, files]
            bam: meta.type == 'bam'
                return [meta, files]    
        }
        .set { ch_reads }

    PREPROCESSING (
        ch_reads.raw
    )
    ch_versions = ch_versions.mix(PREPROCESSING.out.versions)

    //
    // SUBWORKFLOW: GENOME_ASSEMBLY
    // Mix preprocessed reads with already-cleaned reads
    //
    GENOME_ASSEMBLY (
        PREPROCESSING.out.fastp_reads.mix(ch_reads.cleaned)
    )
    ch_versions = ch_versions.mix(GENOME_ASSEMBLY.out.versions)

    //
    // SUBWORKFLOW: Contamination Detection
    // Process genome assemblies from three sources:
    // 1. FASTA files (ch_fasta)
    // 2. BAM files from samplesheet (ch_reads.bam)
    // 3. BAM files from workflow emits (ch_bam)
    //
    CONTAMINATION_DETECTION (
        ch_fasta
            .mix(ch_reads.bam)    // BAMs from samplesheet with type == 'bam'
            .mix(ch_bam),         // BAMs from workflow outputs or other sources
        params.ramdisk_path ?: [],
        params.db_path
    )
    ch_versions = ch_versions.mix(CONTAMINATION_DETECTION.out.versions)
    
    //
    // Collate and save software versions
    //
    def topic_versions = channel.topic("versions")
        .distinct()
        .branch { entry ->
            versions_file: entry instanceof Path
            versions_tuple: true
        }

    def topic_versions_string = topic_versions.versions_tuple
        .map { process, tool, version ->
            [ process[process.lastIndexOf(':')+1..-1], "  ${tool}: ${version}" ]
        }
        .groupTuple(by:0)
        .map { process, tool_versions ->
            tool_versions.unique().sort()
            "${process}:\n${tool_versions.join('\n')}"
        }

    def ch_versions_files = ch_versions.filter { it instanceof Path }

    softwareVersionsToYAML(ch_versions_files.mix(topic_versions.versions_file))
        .mix(topic_versions_string)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name: 'nf_core_'  +  'fsptest_software_'  + 'mqc_'  + 'versions.yml',
            sort: true,
            newLine: true
        ).set { ch_collated_versions }

    //
    // MODULE: MultiQC
    //
    ch_multiqc_config        = channel.fromPath(
        "$projectDir/assets/multiqc_config.yml", checkIfExists: true)
    ch_multiqc_custom_config = params.multiqc_config ?
        channel.fromPath(params.multiqc_config, checkIfExists: true) :
        channel.empty()
    ch_multiqc_logo          = params.multiqc_logo ?
        channel.fromPath(params.multiqc_logo, checkIfExists: true) :
        channel.empty()

    summary_params      = paramsSummaryMap(
        workflow, parameters_schema: "nextflow_schema.json")
    ch_workflow_summary = channel.value(paramsSummaryMultiqc(summary_params))
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_custom_methods_description = params.multiqc_methods_description ?
        file(params.multiqc_methods_description, checkIfExists: true) :
        file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
    ch_methods_description                = channel.value(
        methodsDescriptionText(ch_multiqc_custom_methods_description))

    ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_methods_description.collectFile(
            name: 'methods_description_mqc.yaml',
            sort: true
        )
    )

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList(),
        [],
        []
    )

    emit:
    multiqc_report  = MULTIQC.out.report.toList() // channel: /path/to/multiqc_report.html
    classifications = CONTAMINATION_DETECTION.out.tiara_classifications
    versions        = ch_versions                 // channel: [ path(versions.yml) ]
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/