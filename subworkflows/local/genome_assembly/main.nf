// TODO nf-core: If in doubt look at other nf-core/subworkflows to see how we are doing things! :)
//               https://github.com/nf-core/modules/tree/master/subworkflows
//               You can also ask for help via your pull request or on the #subworkflows channel on the nf-core Slack workspace:
//               https://nf-co.re/join
// TODO nf-core: A subworkflow SHOULD import at least two modules

include { SEQKIT_STATS                         } from '../../../modules/nf-core/seqkit/stats/main'
// include { SEQKIT_STATS as SEQKIT_STATS_MERGED  } from '../../../modules/nf-core/seqkit/stats/main' I can't work on this if the preprocessing subworkflow is not updated to ouput merged reads.
include { FASTK_FASTK                          } from '../../../modules/nf-core/fastk/fastk/main'
include { SPADES                               } from '../../../modules/nf-core/spades/main'
include { MEGAHIT                              } from '../../../modules/nf-core/megahit/main'
include { MINIA                                } from '../../../modules/nf-core/minia/main'
include { RENAME_ASSEMBLIES                    } from '../../../modules/local/rename_assemblies/main' 
include { BUSCO_BUSCO                          } from '../../../modules/nf-core/busco/busco/main'
include { BUSCO_BUSCO as BUSCO_SPECIFIC        } from '../../../modules/nf-core/busco/busco/main'
include { MERQURYFK_MERQURYFK                  } from '../../../modules/nf-core/merquryfk/merquryfk/main'
include { QUAST                                } from '../../../modules/nf-core/quast/main'

workflow GENOME_ASSEMBLY {

    take:
    // TODO nf-core: edit input (take) channels
    ch_fastp_reads // channel: [ val(meta), path(reads) ]

    main:
    // TODO nf-core: substitute modules here for the modules of your subworkflow

    SEQKIT_STATS ( ch_fastp_reads )
//    SEQKIT_STATS_MERGED
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

    // input channel for renaming the assemblies. I need to change the meta.id to include the assembler and avoid conflicts in the output names.
    def ch_draft_assemblies_input = SPADES.out.scaffolds.map { meta, scaffolds -> 
        // add assembler name to meta.id to ensure unique output names
        def assembler = 'spades'
        def new_meta = meta + [assembly_id: "${meta.id}_${assembler}", assembler: 'spades', id: meta.id]
        return [ new_meta, scaffolds, "${meta.id}_${assembler}.fa" ]
    }
    .mix( MEGAHIT.out.contigs.map { meta, contigs -> 
        def assembler = 'megahit'
        def new_meta = meta + [assembly_id: "${meta.id}_${assembler}", assembler: 'megahit', id: meta.id]
        return [ new_meta, contigs, "${meta.id}_${assembler}.fa" ]
    } )
    .mix( MINIA.out.contigs.map { meta, contigs -> 
        def assembler = 'minia'
        def new_meta = meta + [assembly_id: "${meta.id}_${assembler}", assembler: 'minia', id: meta.id]
        return [ new_meta, contigs, "${meta.id}_${assembler}.fa" ]
    } )

    RENAME_ASSEMBLIES ( ch_draft_assemblies_input )

    BUSCO_BUSCO ( RENAME_ASSEMBLIES.out.renamed_assemblies, params.busco_mode, params.busco_lineage, params.busco_lineages_path ?:[], params.busco_config_file ?:[], params.busco_clean_intermediates )
    BUSCO_SPECIFIC ( RENAME_ASSEMBLIES.out.renamed_assemblies, params.busco_mode, params.busco_lineage_specific, params.busco_lineages_path ?:[], params.busco_config_file ?:[], params.busco_clean_intermediates )

    // input channel for the first input required by merquryfk: tuple val(meta) , path(fastk_hist), path(fastk_ktab), path(assembly), path(haplotigs)
    // to obtain this:
    // 1. join fastk hist and ktab in a single list and map to meta.id to be use as key to then join with the assemblies
    def ch_combined_fastk = FASTK_FASTK.out.hist.join(FASTK_FASTK.out.ktab, by: 0).map { meta, hist, ktab -> [ meta.id, hist, ktab ] } 
    // 2. map renamed assemblies to original meta.id (to be used as key to then join with combined fastk results)
    def ch_draft_assemblies_mapped_to_id = RENAME_ASSEMBLIES.out.renamed_assemblies.map { meta, renamed_assembly -> [ meta.id, meta, renamed_assembly ]}
    // 3. join combined fastk with renamed assemblies using meta.id as key 
    def ch_merquryfk_input = ch_combined_fastk.combine( ch_draft_assemblies_mapped_to_id, by: 0 ).map { sample_id, hist, ktab, meta, assembly -> [ meta, hist, ktab, assembly, [] ] }
    
    MERQURYFK_MERQURYFK ( ch_merquryfk_input, [[],[]], [[],[]] ) // no mathernal and pathernal haplotypes for trio mode

    // input channel for quast: I want to run quast once per sample, so I have to group the different assemblies per sample name
    // ch_draft_assemblies_mapped_to_id is: [ sample_id, meta, assembly ]
    // groupTuple(by: 0) groups by position 0 (sample_id)
    // Result: [ sample_id, [meta1, meta2, meta3], [assembly1, assembly2, assembly3] ]
    def ch_quast_input = ch_draft_assemblies_mapped_to_id.groupTuple( by: 0 ).map { sample_id, metas, assemblies -> [ [id: sample_id], assemblies ]}

    QUAST ( ch_quast_input,[[],[]], [[],[]] ) // no reference fasta or gff for quast

    emit:
    // TODO nf-core: edit emitted channels
    seqkit_stats                 = SEQKIT_STATS.out.stats           // channel: [ val(meta), [ bam ] ]
    fastk_ktab                   = FASTK_FASTK.out.ktab             // channel: [ val(meta), path('*.ktab') ]
    fastk_hist                   = FASTK_FASTK.out.hist             // channel: [ val(meta), path('*.hist') ]
    spades_scaffolds             = SPADES.out.scaffolds             // channel: [ val(meta), path('*.scaffolds.fa.gz') ]
    megahit_contigs              = MEGAHIT.out.contigs              // channel: [ val(meta), path('*.contigs.fa.gz') ]
    minia_contigs                = MINIA.out.contigs                // channel: [ val(meta), path('*.contigs.fa') ]
    renamed_assemblies           = RENAME_ASSEMBLIES.out.renamed_assemblies // channel: [ val(meta), path('*.fa.gz') ]
    busco_batch_summary          = BUSCO_BUSCO.out.batch_summary  // channel: [ val(meta), path('*.busco.batch_summary.txt') ]
    busco_short_summaries_txt    = BUSCO_BUSCO.out.short_summaries_txt  // channel: [ val(meta), path('short_summary.*.txt') ]
    merquryfk_completeness_stats = MERQURYFK_MERQURYFK.out.stats // channel: [ val(meta), path('*.completeness.stats') ]
    quast_results                = QUAST.out.results         // channel: [ val(meta), path("${prefix}") ]
}
