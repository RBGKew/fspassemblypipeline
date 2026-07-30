/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { PREPROCESSING               } from '../subworkflows/local/preprocessing/main'
include { GENOME_ASSEMBLY             } from '../subworkflows/local/genome_assembly/main'
include { GENOME_ASSEMBLY_MERGED      } from '../subworkflows/local/genome_assembly_merged/main'
include { SELECT_BEST_ASSEMBLY_AND_QC } from '../subworkflows/local/select_best_assembly_and_qc/main'
include { CONTAMINATION_DETECTION     } from '../subworkflows/local/contamination_detection/main'
include { BLOBTOOLS                   } from '../subworkflows/local/blobtools/main'
include { MULTIQC                     } from '../modules/nf-core/multiqc/main'
include { paramsSummaryMap            } from 'plugin/nf-schema'
include { paramsSummaryMultiqc        } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML      } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText      } from '../subworkflows/local/utils_nfcore_fspassemblypipeline_pipeline'
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow FSPASSEMBLYPIPELINE {

    take:
    ch_samplesheet // channel: samplesheet read in from --input
    multiqc_config
    multiqc_logo
    multiqc_methods_description
    outdir

    main:

    // Validate assembly mode selection
    if (!params.use_paired_reads && !params.use_merged_reads) {
        error "ERROR: At least one assembly mode must be enabled. " +
              "Set --use_paired_reads true or --use_merged_reads true"
    }

    ch_versions = channel.empty()
    ch_multiqc_files = channel.empty()

    //
    // SUBWORKFLOW: Run PREPROCESSING
    //

    ch_samplesheet = ch_samplesheet
        .branch { meta, files ->
            raw: meta.type == 'raw'
            cleaned: meta.type == 'cleaned'
            bam: meta.type == 'bam'
        }

    PREPROCESSING (
        ch_samplesheet.raw
    )

    ch_samplesheet_raw = PREPROCESSING.out.fastp_reads
        .branch { meta, reads ->
            def ranks = (meta.family || meta.order || meta.class || meta.phylum)
            taxonomy: ranks // Run assembly if taxonomy information is given
            no_taxonomy: !ranks // Do not run assembly if no taxonomy information is given
        }

    // Initialise assembly output channels as empty
    def ch_draft_assemblies_paired = channel.empty()
    def ch_draft_assemblies_merged = channel.empty()

    // Conditional: Run paired-end assembly workflow
    if (params.use_paired_reads) {
        GENOME_ASSEMBLY(
            ch_samplesheet_raw.taxonomy.mix(ch_samplesheet.cleaned)
        )
        ch_draft_assemblies_paired = GENOME_ASSEMBLY.out.draft_assemblies_paired
    }

    // Conditional: Run merged reads assembly workflow
    if (params.use_merged_reads) {
        GENOME_ASSEMBLY_MERGED(
            PREPROCESSING.out.fastp_reads_merged,
            PREPROCESSING.out.fastp_reads_unmerged,
            PREPROCESSING.out.fastp_reads
        )
        ch_draft_assemblies_merged = GENOME_ASSEMBLY_MERGED.out.draft_assemblies_merged
    }

    // To avoid triggering SELECT_BEST_ASSEMBLY_AND_QC when type is 'raw' and no taxonomy is provided
    // Combine draft assemblies from both assembly workflows
    ch_assemblies_combined = ch_draft_assemblies_paired.mix(ch_draft_assemblies_merged)
    // Only pass reads for samples that actually have assemblies
    ch_reads_for_select = PREPROCESSING.out.fastp_reads
        .mix(ch_samplesheet.cleaned)
        .map { meta, reads -> [meta.id, meta, reads] }
        .join(
            ch_assemblies_combined
                .map { meta, _assembly, _name -> [meta.id] }
                .unique()
        )
        .map { _id, meta, reads -> [meta, reads] }

    SELECT_BEST_ASSEMBLY_AND_QC(
        ch_draft_assemblies_paired,
        ch_draft_assemblies_merged,
        ch_reads_for_select
    )

    ch_contamination_detection_input = SELECT_BEST_ASSEMBLY_AND_QC.out.best_assembly_polished
    .join(SELECT_BEST_ASSEMBLY_AND_QC.out.bwamem2_best_assembly_bam)
    .join(SELECT_BEST_ASSEMBLY_AND_QC.out.busco_best_assembly_full_table_specific)
    .map { meta, fa, bam, full_table -> [meta, [fa, bam, full_table]] }
    .filter { meta, files -> meta.taxid }  // only run if taxid is set in samplesheet

    CONTAMINATION_DETECTION(
    ch_contamination_detection_input.mix(ch_samplesheet.bam)     // Channel: [meta, fasta, bam, busco]
    )

    BLOBTOOLS(
    ch_contamination_detection_input.mix(ch_samplesheet.bam),
    CONTAMINATION_DETECTION.out.blobtools_taxonomy)

    // ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]})

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

    ch_collated_versions = softwareVersionsToYAML(ch_versions_files.mix(topic_versions.versions_file))
        .mix(topic_versions_string)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name: 'nf_core_'  +  'fsptest_software_'  + 'mqc_'  + 'versions.yml',
            sort: true,
            newLine: true
        )

    //
    // MODULE: MultiQC
    //
    def mqc_config        = file("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
    def mqc_custom_config = params.multiqc_config ? file(params.multiqc_config, checkIfExists: true) : []
    def mqc_logo          = params.multiqc_logo ? file(params.multiqc_logo, checkIfExists: true) : []

    summary_params      = paramsSummaryMap(
        workflow, parameters_schema: "nextflow_schema.json")
    ch_workflow_summary = channel.value(paramsSummaryMultiqc(summary_params))
        .collectFile(name: 'workflow_summary_mqc.yaml')

    ch_multiqc_custom_methods_description = params.multiqc_methods_description ?
        file(params.multiqc_methods_description, checkIfExists: true) :
        file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
    ch_methods_description = channel.value(
        methodsDescriptionText(ch_multiqc_custom_methods_description))
        .collectFile(name: 'methods_description_mqc.yaml')

    // Add global context files (summary, versions, methods) with empty meta
    ch_multiqc_files = ch_multiqc_files
        .mix(ch_workflow_summary.map   { file -> [[:], file] })
        .mix(ch_collated_versions.map  { file -> [[:], file] })
        .mix(ch_methods_description.map{ file -> [[:], file] })

    // Build the MultiQC config list (default + optional custom)
    def mqc_config_files = mqc_custom_config ? [mqc_config, mqc_custom_config] : [mqc_config]

    // Normal merged MultiQC: collect all files into one report
    ch_all_mqc_files = ch_multiqc_files
        .map { _meta, file -> file }
        .collect()

    ch_multiqc_input = ch_all_mqc_files
        .map { files ->
            [
                [id: 'multiqc_report'],
                files,
                mqc_config_files,
                mqc_logo,
                [],  // replace_names
                []   // sample_names
            ]
        }

    MULTIQC(ch_multiqc_input)

    emit:
    multiqc_report = MULTIQC.out.report.map { meta, report -> report }
    versions       = ch_versions



}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
