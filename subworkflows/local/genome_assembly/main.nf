// TODO nf-core: If in doubt look at other nf-core/subworkflows to see how we are doing things! :)
//               https://github.com/nf-core/modules/tree/master/subworkflows
//               You can also ask for help via your pull request or on the #subworkflows channel on the nf-core Slack workspace:
//               https://nf-co.re/join
// TODO nf-core: A subworkflow SHOULD import at least two modules

include { SEQKIT_STATS        } from '../../../modules/nf-core/seqkit/stats/main'
include { FASTK_FASTK         } from '../../../modules/nf-core/fastk/fastk/main'
include { SPADES              } from '../../../modules/nf-core/spades/main'
include { MEGAHIT             } from '../../../modules/nf-core/megahit/main'
include { MINIA               } from '../../../modules/nf-core/minia/main'
include { BUSCO_BUSCO         } from '../../../modules/nf-core/busco/busco/main'
// include { QUAST_QUAST         } from '../../../modules/nf-core/quast/main'
// include { MERQURYFK_MERQURYFK } from '../../../modules/nf-core/merquryfk/merquryfk/main'

workflow GENOME_ASSEMBLY {

    take:
    // TODO nf-core: edit input (take) channels
    ch_fastp_reads // channel: [ val(meta), path(reads) ]

    main:
    // TODO nf-core: substitute modules here for the modules of your subworkflow

    SEQKIT_STATS ( ch_fastp_reads )
    FASTK_FASTK  ( ch_fastp_reads )

    // Spades needs a tuple with 4 elements as inputs, so we need to map the channel to add empty lists for the other 2 inputs (see PREPROCESSING subworkflow for example)
    // SPADES: [ meta, illumina, pacbio, nanopore ]
    ch_input_reads_spades = ch_fastp_reads.map { meta, reads -> [ meta, reads, [], [] ] }

    SPADES       ( ch_input_reads_spades,
    [],
    []
    )

    // Megahit needs a tuple with 3 elements as input. I can't use ch_fastp_reads directly because R1 and R2 paths there are in a single list element. So I need to map the channel to split R1 and R2 into separate list elements.
    // MEGAHIT: [ meta, reads1, reads2 ]
    ch_input_reads_megahit = ch_fastp_reads.map { meta, reads -> [ meta, reads[0], reads[1] ] }

    MEGAHIT      ( ch_input_reads_megahit )

    MINIA        ( ch_fastp_reads )

    // input channel for BUSCO
    def ch_busco_input = SPADES.out.scaffolds.map { meta, scaffolds -> 
        // add assembler name to meta.id to ensure unique output names
        def new_meta = meta + [id: "${meta.id}_spades", assembler: 'spades']
        return [ new_meta, scaffolds ]
    }
    .mix( MEGAHIT.out.contigs.map { meta, contigs -> 
        def new_meta = meta + [id: "${meta.id}_megahit", assembler: 'megahit']
        return [ new_meta, contigs ]
    } )
    .mix( MINIA.out.contigs.map { meta, contigs -> 
        def new_meta = meta + [id: "${meta.id}_minia", assembler: 'minia']
        return [ new_meta, contigs ]
    } )

    BUSCO_BUSCO ( ch_busco_input, params.busco_mode, params.busco_lineage, params.busco_lineages_path ?:[], params.busco_config_file ?:[], params.busco_clean_intermediates )

     // For Merqury we need to provide a list, so rather than splitting the channel, this time we need to put things together into a list element.
     // MERQURYFK: [ meta, fastk_hist, fastk_ktab, assembly, haplotigs ]

//    MERQURYFK_MERQURYFK ( FASTK_FASTK.out.hist, FASTK_FASTK.out.ktab, ch_busco_input )

    emit:
    // TODO nf-core: edit emitted channels
    seqkit_stats      = SEQKIT_STATS.out.stats           // channel: [ val(meta), [ bam ] ]
    fastk_ktab        = FASTK_FASTK.out.ktab             // channel: [ val(meta), path('*.ktab') ]
    fastk_hist        = FASTK_FASTK.out.hist             // channel: [ val(meta), path('*.hist') ]
    spades_scaffolds  = SPADES.out.scaffolds             // channel: [ val(meta), path('*.scaffolds.fa.gz') ]
    megahit_contigs   = MEGAHIT.out.contigs              // channel: [ val(meta), path('*.contigs.fa.gz') ]
    minia_contigs     = MINIA.out.contigs                // channel: [ val(meta), path('*.contigs.fa') ]
    busco_batch_summary = BUSCO_BUSCO.out.batch_summary  // channel: [ val(meta), path('*.busco.batch_summary.txt') ]
    busco_short_summaries_txt = BUSCO_BUSCO.out.short_summaries_txt  // channel: [ val(meta), path('short_summary.*.txt') ]
//    merquryfk_completeness_stats = MERQURYFK_MERQURYFK.out.stats // channel: [ val(meta), path('*.completeness.stats') ]
}
