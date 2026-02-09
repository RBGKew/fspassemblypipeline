// TODO nf-core: If in doubt look at other nf-core/subworkflows to see how we are doing things! :)
//               https://github.com/nf-core/modules/tree/master/subworkflows
//               You can also ask for help via your pull request or on the #subworkflows channel on the nf-core Slack workspace:
//               https://nf-co.re/join
// TODO nf-core: A subworkflow SHOULD import at least two modules

include { SEQKIT_STATS      } from '../../../modules/nf-core/seqkit/stats/main'
// include { SAMTOOLS_INDEX     } from '../../../modules/nf-core/samtools/index/main'

workflow GENOME_ASSEMBLY {

    take:
    // TODO nf-core: edit input (take) channels
    ch_fastp_reads // channel: [ val(meta), path(reads) ]

    main:
    // TODO nf-core: substitute modules here for the modules of your subworkflow

    SEQKIT_STATS ( ch_fastp_reads )

//    SAMTOOLS_INDEX ( SAMTOOLS_SORT.out.bam )

    emit:
    // TODO nf-core: edit emitted channels
    seqkit_stats      = SEQKIT_STATS.out.stats           // channel: [ val(meta), [ bam ] ]
}
