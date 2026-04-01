#!/usr/bin/env nextflow

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Local debugging entrypoint for PREPROCESSING subworkflow only
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { PREPROCESSING  } from './subworkflows/local/preprocessing/main'
include { samplesheetToList } from 'plugin/nf-schema'

workflow {

    main:
    channel
        .fromList(samplesheetToList(params.input, "${projectDir}/assets/schema_input.json"))
        .map {
            meta, file_1, file_2, fasta ->
                def sample_type_key = "${meta.id}__${meta.type}"
                if ( meta.type == 'bam') {
                    return [ sample_type_key, meta + [ single_end:false ], [ fasta, file_1 ] ]
                } else if (!file_2) {
                    return [ sample_type_key, meta + [ single_end:true ], [ file_1 ] ]
                } else {
                    return [ sample_type_key, meta + [ single_end:false ], [ file_1, file_2 ] ]
                }
        }
        .groupTuple()
        .map { samplesheet -> validateInputSamplesheet(samplesheet) }
        .map { meta, fastqs -> [ meta, fastqs.flatten() ] }
        .branch { meta, files ->
            raw: meta.type == 'raw'
            cleaned: meta.type == 'cleaned'
            bam: meta.type == 'bam'
        }
        .set { ch_samplesheet }

    PREPROCESSING(
        ch_samplesheet.raw
    )
}

def validateInputSamplesheet(input) {
    def (metas, fastqs) = input[1..2]

    def endedness_ok = metas.collect { meta -> meta.single_end }.unique().size == 1
    if (!endedness_ok) {
        error("Please check input samplesheet -> Multiple runs of a sample must be of the same datatype i.e. single-end or paired-end: ${metas[0].id}")
    }

    return [ metas[0], fastqs ]
}
