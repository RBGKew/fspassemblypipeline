// TODO nf-core: If in doubt look at other nf-core/subworkflows to see how we are doing things! :)
//               https://github.com/nf-core/modules/tree/master/subworkflows
//               You can also ask for help via your pull request or on the #subworkflows channel on the nf-core Slack workspace:
//               https://nf-co.re/join
// TODO nf-core: A subworkflow SHOULD import at least two modules

include { FALCO      } from '../../../modules/nf-core/falco/main'
include { FASTQC     } from '../../../modules/nf-core/fastqc/main'
include { FASTP      } from '../../../modules/nf-core/fastp/main'

workflow PREPROCESSING {

    take:

    ch_samplesheet // channel: [ val(meta), [ reads ] ]

    main:

    ch_versions = channel.empty()

    FASTQC (
        ch_samplesheet
    )
    ch_versions = ch_versions.mix( FASTQC.out.versions )

    // Prepare channel for FASTP
    ch_samplesheet_fastp = ch_samplesheet.map { meta, reads -> [ meta, reads, [] ] }

    FASTP (
        ch_samplesheet_fastp,
        [],
        [],
        []
    )
    // ch_versions = ch_versions.mix( FASTP.out.versions_fastp )

    FALCO (
        FASTP.out.reads
    )
    ch_versions = ch_versions.mix( FALCO.out.versions )

    emit:
    // Emit module versions
    versions = ch_versions
}