#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
<<<<<<< HEAD
    nf-core/fspassemblypipeline
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/nf-core/fspassemblypipeline
    Website: https://nf-co.re/fspassemblypipeline
    Slack  : https://nfcore.slack.com/channels/fspassemblypipeline
=======
    nf-core/fsptest
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/nf-core/fsptest
    Website: https://nf-co.re/fsptest
    Slack  : https://nfcore.slack.com/channels/fsptest
>>>>>>> NiallG1/fsptest/dev
----------------------------------------------------------------------------------------
*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS / WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

<<<<<<< HEAD
include { FSPASSEMBLYPIPELINE  } from './workflows/fspassemblypipeline'
include { PIPELINE_INITIALISATION } from './subworkflows/local/utils_nfcore_fspassemblypipeline_pipeline'
include { PIPELINE_COMPLETION     } from './subworkflows/local/utils_nfcore_fspassemblypipeline_pipeline'
=======
include { FSPTEST  } from './workflows/fsptest'
include { PIPELINE_INITIALISATION } from './subworkflows/local/utils_nfcore_fsptest_pipeline'
include { PIPELINE_COMPLETION     } from './subworkflows/local/utils_nfcore_fsptest_pipeline'
>>>>>>> NiallG1/fsptest/dev
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    NAMED WORKFLOWS FOR PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// WORKFLOW: Run main analysis pipeline depending on type of input
//
<<<<<<< HEAD
workflow NFCORE_FSPASSEMBLYPIPELINE {

    take:
    samplesheet // channel: samplesheet read in from --input
=======
workflow NFCORE_FSPTEST {

    take:
    samplesheet // channel: samplesheet read in from --input
    tiara_input
>>>>>>> NiallG1/fsptest/dev

    main:

    //
    // WORKFLOW: Run pipeline
    //
<<<<<<< HEAD
    FSPASSEMBLYPIPELINE (
        samplesheet
    )
    emit:
    multiqc_report = FSPASSEMBLYPIPELINE.out.multiqc_report // channel: /path/to/multiqc_report.html
=======
    FSPTEST (
        samplesheet,
        tiara_input
    )
    emit:
    multiqc_report = FSPTEST.out.multiqc_report // channel: /path/to/multiqc_report.html
    classifications = FSPTEST.out.classifications
>>>>>>> NiallG1/fsptest/dev
}
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow {

    main:
    //
    // SUBWORKFLOW: Run initialisation tasks
    //
    PIPELINE_INITIALISATION (
        params.version,
        params.validate_params,
        params.monochrome_logs,
        args,
        params.outdir,
        params.input,
<<<<<<< HEAD
=======
        params.tiara_input,
>>>>>>> NiallG1/fsptest/dev
        params.help,
        params.help_full,
        params.show_hidden
    )

    //
    // WORKFLOW: Run main workflow
    //
<<<<<<< HEAD
    NFCORE_FSPASSEMBLYPIPELINE (
        PIPELINE_INITIALISATION.out.samplesheet
=======
    NFCORE_FSPTEST (
        PIPELINE_INITIALISATION.out.samplesheet,
        PIPELINE_INITIALISATION.out.tiara_input
>>>>>>> NiallG1/fsptest/dev
    )
    //
    // SUBWORKFLOW: Run completion tasks
    //
    PIPELINE_COMPLETION (
        params.email,
        params.email_on_fail,
        params.plaintext_email,
        params.outdir,
        params.monochrome_logs,
        params.hook_url,
<<<<<<< HEAD
        NFCORE_FSPASSEMBLYPIPELINE.out.multiqc_report
=======
        NFCORE_FSPTEST.out.multiqc_report
>>>>>>> NiallG1/fsptest/dev
    )
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
