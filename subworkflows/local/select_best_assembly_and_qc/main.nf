include { FASTK_FASTK as FASTK_ASM                         } from '../../../modules/nf-core/fastk/fastk/main'
include { RENAME_ASSEMBLIES                                } from '../../../modules/local/rename_assemblies/main'
include { BUSCO_BUSCO                                      } from '../../../modules/nf-core/busco/busco/main'
include { BUSCO_BUSCO as BUSCO_SPECIFIC                    } from '../../../modules/nf-core/busco/busco/main'
include { MERQURYFK_MERQURYFK                              } from '../../../modules/nf-core/merquryfk/merquryfk/main'
include { QUAST                                            } from '../../../modules/nf-core/quast/main'
include { SELECTBESTASSEMBLY                               } from '../../../modules/local/selectbestassembly/main'
include { PYPOLCA_RUN                                      } from '../../../modules/nf-core/pypolca/run/main'
include { BUSCO_BUSCO as BUSCO_FINAL                       } from '../../../modules/nf-core/busco/busco/main'
include { BUSCO_BUSCO as BUSCO_SPECIFIC_FINAL              } from '../../../modules/nf-core/busco/busco/main'
include { MERQURYFK_MERQURYFK as MERQURYFK_FINAL           } from '../../../modules/nf-core/merquryfk/merquryfk/main'
include { QUAST as QUAST_FINAL                             } from '../../../modules/nf-core/quast/main'
include { BWAMEM2_INDEX                                    } from '../../../modules/nf-core/bwamem2/index/main'
include { BWAMEM2_MEM                                      } from '../../../modules/nf-core/bwamem2/mem/main'
include { SAMTOOLS_INDEX                                   } from '../../../modules/nf-core/samtools/index/main'
include { SAMTOOLS_FAIDX                                   } from '../../../modules/nf-core/samtools/faidx/main'
include { SAMTOOLS_COVERAGE                                } from '../../../modules/nf-core/samtools/coverage/main'
include { SAMTOOLS_FLAGSTAT                                } from '../../../modules/nf-core/samtools/flagstat/main'
include { COVERAGEVIZ                                      } from '../../../modules/local/coverageviz/main'
include { EXTRACTASSEMBLYMETRICS                            } from '../../../modules/local/extractassemblymetrics/main'

workflow SELECT_BEST_ASSEMBLY_AND_QC {

    take:
    ch_draft_assemblies_paired // channel: draft assemblies from genome_assembly
    ch_draft_assemblies_merged // channel: draft assemblies from genome_assembly_merged
    ch_fastp_reads // needed to run fastk. We run fastk only with R1R2 reads, even for the merged assembly, so that we can see errors introduced by the merging step.

    main:


    FASTK_ASM  ( ch_fastp_reads )

// =================== Draft assemblies QC =======================

    // input to rename_assemblies all the assemblies coming from both assembly workflows

    ch_draft_assemblies_input = ch_draft_assemblies_paired
        .mix(ch_draft_assemblies_merged)

    RENAME_ASSEMBLIES ( ch_draft_assemblies_input )

    BUSCO_BUSCO ( RENAME_ASSEMBLIES.out.renamed_assemblies, params.busco_mode, params.busco_lineage, params.busco_lineages_path ?:[], params.busco_config_file ?:[], params.busco_clean_intermediates )

    // Map samples to their best BUSCO database based on taxonomy
    // BUSCO_SPECIFIC runs BUSCO with a specific lineage for each sample. The lineage is determined based on the sample metadata and the list of available lineages.
    // we attempt to find the most specific available BUSCO database for each sample by walking down the taxonomic hierarchy (family → order → class → phylum).
    // If no matching database is found at any level, we fall back to the lineage specified by params.busco_lineage.

    // Load BUSCO lineages once
    def busco_lineages_file = file(params.lineages_list_file)

    def busco_lineages = busco_lineages_file
        .readLines() as Set

    def ext = "_${params.busco_db_extension}"
    def fallback = params.busco_lineage

    // Add lineage to each sample and split into synchronized channels
    ch_assemblies_with_lineage = RENAME_ASSEMBLIES.out.renamed_assemblies
        .map { meta, assembly ->
            // Find best matching lineage for this sample
            def lineage = ['family', 'order', 'class', 'phylum']
                .findResult { level ->
                    def taxon = meta[level]
                    if (!taxon) {
                        log.info "  ${level}: not set"
                        return null
                    }
                    def candidate = "${taxon.toLowerCase()}${ext}".toString()
                    def found = candidate in busco_lineages
                    found ? candidate : null
                } ?: fallback

            // Return tuple with lineage for splitting
            tuple(meta, assembly, lineage)
        }
        .multiMap { meta, assembly, lineage ->
            assemblies: tuple(meta, assembly)
            lineages: lineage
        }

    // Specific BUSCO with per-sample lineage
    BUSCO_SPECIFIC(
        ch_assemblies_with_lineage.assemblies,
        params.busco_mode,
        ch_assemblies_with_lineage.lineages,
        params.busco_lineages_path ?: [],
        params.busco_config_file ?: [],
        params.busco_clean_intermediates
    )



    // input channel for the first input required by merquryfk: tuple val(meta) , path(fastk_hist), path(fastk_ktab), path(assembly), path(haplotigs)
    // to obtain this:
    // 1. join fastk hist and ktab in a single list and map to meta.id to be use as key to then join with the assemblies
    def ch_combined_fastk = FASTK_ASM.out.hist.join(FASTK_ASM.out.ktab, by: 0).map { meta, hist, ktab -> [ meta.id, hist, ktab ] }
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

// =================== Select best assembly =======================

    // Group BUSCO summaries per sample
    def ch_busco_per_sample = BUSCO_SPECIFIC.out.short_summaries_txt
        .map { meta, summary -> [ meta.id, summary ] }
        .groupTuple(by: 0)
    // [ sample_id, [summary1, summary2, ...] ]

    // QUAST.out.tsv already emits [ [id: sample_id], tsv ]
    // map to plain id for joining
    def ch_quast_tsv = QUAST.out.tsv
        .map { meta, tsv -> [ meta.id, tsv ] }
    // [ sample_id, tsv ]

    // Group assemblies per sample
    def ch_assemblies_per_sample = RENAME_ASSEMBLIES.out.renamed_assemblies
        .map { meta, asm -> [ meta.id, asm ] }
        .groupTuple(by: 0)
    // [ sample_id, [asm1, asm2, ...] ]


    def ch_select_best_assembly_input = ch_busco_per_sample
        .join(ch_quast_tsv,          by: 0)
        .join(ch_assemblies_per_sample, by: 0)
        .map { sample_id, buscos, tsv, asms ->
            [ [id: sample_id], buscos, tsv, asms ]
        }

    // Select best assembly based on BUSCO and QUAST metrics. The output is a channel with the best assembly per sample, along with its metadata and metrics.

    SELECTBESTASSEMBLY ( ch_select_best_assembly_input )

    // re-add meta information to the selected best assembly channel for downstream use
    def ch_selected_best_assembly = SELECTBESTASSEMBLY.out.best_assembly
        .join(SELECTBESTASSEMBLY.out.best_assembly_meta, by: 0)
        .map { meta, fa, meta_file ->
            def fields = meta_file.text.readLines()
                .collectEntries { line -> line.split('=') as List }
            [ meta + fields, fa ]
        }
    // emits: [[id:'sample_id', reads_type:'R1R2', kmer_strategy:'kmergenie', assembler:'spades'], fa]

// =================== Best assembly correction with pypolca =======================

    // Pypolca needs reads and assembly paired by meta.id. However pypolca expects two tuples as input with meta and path
    // so we need to sync them by meta.id and then split them again into the two tuples expected by pypolca

    def ch_fastp_reads_by_id = ch_fastp_reads.map { meta, reads -> [meta.id, meta, reads] }
    def ch_selected_best_assembly_by_id = ch_selected_best_assembly.map { meta, fa -> [meta.id, meta, fa] }

    // Also make sure that the meta of the assembly is carried forward.
    // Since I need to add some steps, I'll carry forward the full meta (assembly+reads)

    // Join reads and assembly by meta.id
    def ch_reads_asm_joined = ch_selected_best_assembly_by_id.join(ch_fastp_reads_by_id, by: 0)
    .map { id, meta_fa, fa, meta_reads, reads -> [meta_fa + meta_reads, fa, reads] }

    // Split into the two tuples expected by pypolca: one with meta and reads, the other with meta and assembly
    def ch_pypolca_reads = ch_reads_asm_joined.map { meta_merged, fa, reads -> [meta_merged, reads] }
    def ch_pypolca_asm   = ch_reads_asm_joined.map { meta_merged, fa, reads -> [meta_merged, fa] }

    PYPOLCA_RUN (
        ch_pypolca_reads,
        ch_pypolca_asm
        )

// =================== Best assembly QC =======================

    // Run BUSCO on the polished best assembly using the general lineage specified by params.busco_lineage and then run BUSCO with the specific lineage determined for each sample based on its metadata.

    BUSCO_FINAL ( PYPOLCA_RUN.out.polished, params.busco_mode, params.busco_lineage, params.busco_lineages_path ?:[], params.busco_config_file ?:[], params.busco_clean_intermediates )

    // Add lineage to each sample and split into synchronized channels
    ch_best_assembly_with_lineage = PYPOLCA_RUN.out.polished
        .map { meta, assembly ->
            // Find best matching lineage for this sample
            def lineage = ['family', 'order', 'class', 'phylum']
                .findResult { level ->
                    def taxon = meta[level]
                    if (!taxon) {
                        log.info "  ${level}: not set"
                        return null
                    }
                    def candidate = "${taxon.toLowerCase()}${ext}".toString()
                    def found = candidate in busco_lineages
                    found ? candidate : null
                } ?: fallback

            // Return tuple with lineage for splitting
            tuple(meta, assembly, lineage)
        }
        .multiMap { meta, assembly, lineage ->
            assembly: tuple(meta, assembly)
            lineages: lineage
        }

    BUSCO_SPECIFIC_FINAL(
        ch_best_assembly_with_lineage.assembly,
        params.busco_mode,
        ch_best_assembly_with_lineage.lineages,
        params.busco_lineages_path ?: [],
        params.busco_config_file ?: [],
        params.busco_clean_intermediates
    )

    // Run merquryfk on the polished best assembly using the fastk hist and ktab from the reads.
    //The input channel for merquryfk is a tuple with meta, hist, ktab, polished_assembly, haplotigs (empty list since we don't have haplotigs).

    def ch_polished_assembly_mapped_to_id = PYPOLCA_RUN.out.polished.map { meta, polished_assembly -> [ meta.id, meta, polished_assembly ] }

    def ch_merquryfk_final_input = ch_combined_fastk.combine( ch_polished_assembly_mapped_to_id, by: 0 )
        .map { sample_id, hist, ktab, meta, polished_assembly -> [ meta, hist, ktab, polished_assembly, [] ] }

    MERQURYFK_FINAL ( ch_merquryfk_final_input, [[],[]], [[],[]] )

    // Run quast on the polished best assembly. This time we only have one assembly per sample as input, rather than a list of assemblies.

    QUAST_FINAL ( PYPOLCA_RUN.out.polished,[[],[]], [[],[]] ) // no reference fasta or gff for quast

    // Run bwamem2 to map the reads to the polished best assembly. This is to perform a coverage analysis.

    BWAMEM2_INDEX ( PYPOLCA_RUN.out.polished )

    // map index to meta.id for joining with reads and assembly
    ch_index_mapped_to_id = BWAMEM2_INDEX.out.index.map { meta, index -> [ meta.id, meta, index ] }

    // Join reads, polished assembly and index by meta.id to make sure they are sync and use the meta from the assembly
    def ch_bwamem2_joined = ch_fastp_reads_by_id
        .join(ch_polished_assembly_mapped_to_id, by: 0)
        .join(ch_index_mapped_to_id, by: 0)
        .map { id, meta_reads, reads, meta_asm, assembly, meta_index, index ->
            [meta_asm, reads, assembly, index]
        }

    // Split into the three tuples expected by bwamem2: one with meta and reads, one with meta and index, one with meta and assembly
    def ch_bwamem2_reads = ch_bwamem2_joined.map { meta, reads, assembly, index -> [meta, reads] }
    def ch_bwamem2_index = ch_bwamem2_joined.map { meta, reads, assembly, index -> [meta, index] }
    def ch_bwamem2_asm   = ch_bwamem2_joined.map { meta, reads, assembly, index -> [meta, assembly] }

    BWAMEM2_MEM ( ch_bwamem2_reads,
        ch_bwamem2_index,
        ch_bwamem2_asm,
        true // sort_bam set to true will view and sort the sam file, result: sorted bam file
        )

    // Use samtools to extract coverage and flagstat metrics from the bam file produced by bwamem2.

    SAMTOOLS_INDEX ( BWAMEM2_MEM.out.bam )

    def ch_faidx_input = PYPOLCA_RUN.out.polished.map { meta, polished_assembly -> [meta, polished_assembly, []] }

    SAMTOOLS_FAIDX (
        ch_faidx_input,
        false // do not produce a file containing chromosome sizes.
     )

    // Samtools coverage needs:
    // tuple 1: meta, sorted_bam, index_sorted_bam
    // tuple 2: meta, polished_assembly, fai_polished_assembly

    def ch_coverage_input_bam = BWAMEM2_MEM.out.bam.join(SAMTOOLS_INDEX.out.index, by: 0)
        .map { meta, bam, index -> [meta, bam, index] }

    def ch_coverage_input_fasta = PYPOLCA_RUN.out.polished.join(SAMTOOLS_FAIDX.out.fai, by: 0)
        .map { meta, polished_assembly, fai -> [meta, polished_assembly, fai] }

    // the meta of these two channels is exactly the same so we can use it directly in samtools coverage without joining again

    SAMTOOLS_COVERAGE (
        ch_coverage_input_bam,
        ch_coverage_input_fasta
    )

    // Samtools flagstats only needs one tuple containig meta, bam, index.
    // So we can use the bam channel used for samtools coverage

    SAMTOOLS_FLAGSTAT ( ch_coverage_input_bam )

    // Run coverage_visualisation.py to visualise the coverage and flagstat metrics produced by samtools.
    // Coverageviz needs a tuple with meta, samtools coverage output (txt) and samtools flagstat output (txt)
    def ch_coverageviz_input = SAMTOOLS_COVERAGE.out.coverage.join(SAMTOOLS_FLAGSTAT.out.flagstat, by: 0)
        .map { meta, coverage, flagstat -> [meta, coverage, flagstat] }

    COVERAGEVIZ ( ch_coverageviz_input )

    // Join all per-sample final QC outputs by meta.id
    def ch_metrics_input = BUSCO_FINAL.out.short_summaries_txt
        .join(BUSCO_SPECIFIC_FINAL.out.short_summaries_txt, by: 0)
        .join(QUAST_FINAL.out.results, by: 0)
        .join(MERQURYFK_FINAL.out.qv, by: 0)
        .join(MERQURYFK_FINAL.out.stats, by: 0)
        .join(COVERAGEVIZ.out.summary, by: 0)

    EXTRACTASSEMBLYMETRICS( ch_metrics_input )

    // Merge all per-sample rows into one CSV (collectFile handles the header)
    def ch_all_metrics = EXTRACTASSEMBLYMETRICS.out.metrics_row.collectFile(
        name: 'assembly_metrics.tsv',
        keepHeader: true,
        skip: 1,  // skip header in all but first file
        sort: true, // sort by sample id
        storeDir: "${params.outdir}/genome_assembly"
    )

    emit:
    // Fastk outputs (needed as input for merquryfk)
    fastk_ktab                                           = FASTK_ASM.out.ktab
    fastk_hist                                           = FASTK_ASM.out.hist

    // Draft assemblies QC
    renamed_assemblies                                   = RENAME_ASSEMBLIES.out.renamed_assemblies
    busco_batch_summary                                  = BUSCO_BUSCO.out.batch_summary
    busco_short_summaries_txt                            = BUSCO_BUSCO.out.short_summaries_txt
    busco_full_table                                     = BUSCO_BUSCO.out.full_table
    busco_batch_summary_specific                         = BUSCO_SPECIFIC.out.batch_summary
    busco_short_summaries_txt_specific                   = BUSCO_SPECIFIC.out.short_summaries_txt
    busco_full_table_specific                            = BUSCO_SPECIFIC.out.full_table
    merquryfk_completeness_stats                         = MERQURYFK_MERQURYFK.out.stats
    quast_results                                        = QUAST.out.results

    // Selected best assembly and its QC
    best_assembly_label                                  = SELECTBESTASSEMBLY.out.best_assembly_label
    best_assembly_meta                                   = SELECTBESTASSEMBLY.out.best_assembly_meta
    best_assembly_busco_scores                           = SELECTBESTASSEMBLY.out.busco_scores
    best_assembly_aun_scores                             = SELECTBESTASSEMBLY.out.aun_scores
    best_assembly_polished                               = PYPOLCA_RUN.out.polished
    busco_best_assembly_batch_summary                    = BUSCO_FINAL.out.batch_summary
    busco_best_assembly_short_summaries_txt              = BUSCO_FINAL.out.short_summaries_txt
    busco_best_assembly_full_table                       = BUSCO_FINAL.out.full_table
    busco_best_assembly_batch_summary_specific           = BUSCO_SPECIFIC_FINAL.out.batch_summary
    busco_best_assembly_short_summaries_txt_specific     = BUSCO_SPECIFIC_FINAL.out.short_summaries_txt
    busco_best_assembly_full_table_specific              = BUSCO_SPECIFIC_FINAL.out.full_table
    quast_best_assembly_results                          = QUAST_FINAL.out.results
    merquryfk_best_assembly_completeness_stats           = MERQURYFK_FINAL.out.stats
    bwamem2_best_assembly_bam                            = BWAMEM2_MEM.out.bam
    samtools_best_assembly_coverage                      = SAMTOOLS_COVERAGE.out.coverage
    samtools_best_assembly_flagstat                      = SAMTOOLS_FLAGSTAT.out.flagstat
    coverageviz_best_assembly_coverage_plot              = COVERAGEVIZ.out.png
    coverageviz_best_assembly_coverage_summary           = COVERAGEVIZ.out.summary
    assembly_metrics                                     = ch_all_metrics
}
