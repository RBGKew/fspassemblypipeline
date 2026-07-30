include { SEQKIT_STATS                                     } from '../../../modules/nf-core/seqkit/stats/main'
include { GETSEQKITK                                       } from '../../../modules/local/getseqkitk/main'
include { KMERGENIE                                        } from '../../../modules/nf-core/kmergenie/main'
include { GETKMERGENIEK                                    } from '../../../modules/local/getkmergeniek/main'
include { SPADES_MERGED as SPADES_MANUAL                   } from '../../../modules/local/spades_merged/main'
include { SPADES_MERGED as SPADES_KMERGENIE                } from '../../../modules/local/spades_merged/main'
include { SPADES_MERGED as SPADES_READS_LENGTH             } from '../../../modules/local/spades_merged/main'
include { MEGAHIT as MEGAHIT_MANUAL                        } from '../../../modules/nf-core/megahit/main'
include { MEGAHIT as MEGAHIT_KMERGENIE                     } from '../../../modules/nf-core/megahit/main'
include { MEGAHIT as MEGAHIT_READS_LENGTH                  } from '../../../modules/nf-core/megahit/main'
include { MINIA as MINIA_MANUAL                            } from '../../../modules/nf-core/minia/main'
include { MINIA as MINIA_KMERGENIE                         } from '../../../modules/nf-core/minia/main'
include { MINIA as MINIA_READS_LENGTH                      } from '../../../modules/nf-core/minia/main'
include { ABYSS_ABYSSPE as ABYSS_MANUAL                    } from '../../../modules/nf-core/abyss/abysspe/main'
include { ABYSS_ABYSSPE as ABYSS_KMERGENIE                 } from '../../../modules/nf-core/abyss/abysspe/main'
include { ABYSS_ABYSSPE as ABYSS_READS_LENGTH              } from '../../../modules/nf-core/abyss/abysspe/main'
include { SPARSEASSEMBLER as SPARSEASSEMBLER_MANUAL        } from '../../../modules/local/sparseassembler/main'
include { SPARSEASSEMBLER as SPARSEASSEMBLER_KMERGENIE     } from '../../../modules/local/sparseassembler/main'
include { SPARSEASSEMBLER as SPARSEASSEMBLER_READS_LENGTH  } from '../../../modules/local/sparseassembler/main'
include { MASURCA as MASURCA_MANUAL                        } from '../../../modules/local/masurca/main'
include { MASURCA as MASURCA_KMERGENIE                     } from '../../../modules/local/masurca/main'
include { MASURCA as MASURCA_READS_LENGTH                  } from '../../../modules/local/masurca/main'
include {createAssemblyMeta                                } from '../utils_nfcore_fspassemblypipeline_pipeline/main'

workflow GENOME_ASSEMBLY_MERGED {

    take:
    ch_fastp_reads_merged // channel: [ val(meta), path(reads) ]
    ch_fastp_reads_unmerged // channel: [ val(meta), path(reads) ] - unmerged reads after merge attempt - needed for spades
    ch_fastp_reads // channel: [ val(meta), path(reads) ] - original paired reads, needed by abyss

    main:

    // Add reads_type to the meta
    ch_merged_reads = ch_fastp_reads_merged.map { meta, reads ->
        def new_meta = meta + [ reads_type: 'merged', single_end: true ]
        tuple(new_meta, reads)
    }

// ==================== K-mer strategies for genome assembly =======================

    def ch_manual_strategy = Channel.empty()
    def ch_kmergenie_strategy = Channel.empty()
    def ch_reads_length_strategy = Channel.empty()

    // Channel 1: Manual strategy (uses config values)
    if (!params.skip_manual_strategy) {
        ch_manual_strategy = ch_merged_reads
            .map { meta, reads ->
                [meta + [kmer_strategy: 'manual'], reads]
            }
    }
    // Channel 2: KmerGenie strategy (adds predicted kmer to default list)
    // KMERGENIE and GETKMERGENIEK only run if skip_kmergenie_strategy is false.
    // The channel consequently is only populated if skip_kmergenie_strategy is false.

    // Initialise channels for outputs as empty to avoid use of conditionals in emit section.
    def ch_kmergenie_html = channel.empty()
    def ch_getkmergeniek_k = channel.empty()

    if (!params.skip_kmergenie_strategy) {
        KMERGENIE(ch_merged_reads)
        GETKMERGENIEK(KMERGENIE.out.html)

        // Capture output for emits (avoids using conditionals in emits section)
        ch_kmergenie_html = KMERGENIE.out.html
        ch_getkmergeniek_k = GETKMERGENIEK.out.kmer_txt

        ch_kmergenie_strategy = ch_merged_reads
            .map { meta, reads -> [meta.id, meta, reads] }
            .join(GETKMERGENIEK.out.kmer_txt.map { meta, kmer_file ->
                [meta.id, kmer_file.text.trim() as Integer]
            })
            .map { id, meta, reads, kmergenie_kmer ->
                // Build k-mer list: add predicted kmer to defaults if valid
                def default_kmers_spades = [21, 33, 55, 77]  // spades recommended defaults
                def max_kmer_spades = 127  // spades max kmer limit

                def default_kmers_megahit = [21, 29, 39, 59, 79, 99, 119, 141] // megahit recommended defaults
                def max_kmer_megahit = 141 // megahit max kmer limit

                // Validate predicted kmer for spades: must be odd, in range [15, 127], and not already in list
                def is_valid_spades = (kmergenie_kmer >= 15 &&
                                kmergenie_kmer <= max_kmer_spades &&
                                kmergenie_kmer % 2 == 1 &&
                                !(kmergenie_kmer in default_kmers_spades))

                // Validate predicted kmer for megahit: must be odd, in range [15, 141], and not already in list
                def is_valid_megahit = (kmergenie_kmer >= 15 &&
                                kmergenie_kmer <= max_kmer_megahit &&
                                kmergenie_kmer % 2 == 1 &&
                                !(kmergenie_kmer in default_kmers_megahit))

                // Build final k-mer lists for spades
                def kmer_list_spades = is_valid_spades ?
                    (default_kmers_spades + [kmergenie_kmer]).sort() :
                    default_kmers_spades

                // Build final k-mer lists for megahit
                def kmer_list_megahit = is_valid_megahit ?
                    (default_kmers_megahit + [kmergenie_kmer]).sort() :
                    default_kmers_megahit

                // For assemblers that only take a single k-mer, use the predicted k-mer if it's valid (between 15 and 127), or fall back to a fixed value (25)
                def single_kmer = (kmergenie_kmer >= 15 && kmergenie_kmer <= 127 && kmergenie_kmer % 2 == 1) ? kmergenie_kmer : 25

                // Enrich metadata with k-mer strategy and lists
                def enriched_meta = meta + [
                    kmer_strategy: 'kmergenie',
                    predicted_kmer: kmergenie_kmer,
                    kmer_list_spades: kmer_list_spades.join(','),
                    kmer_list_megahit: kmer_list_megahit.join(','),
                    single_kmer: single_kmer
                ]
                [enriched_meta, reads]
        }
    }

    // Channel 3: reads_length strategy (adds kmer calculated from median reads length to default list)
    // SEQKIT_STATS, GETSEQKITK only run if skip_reads_length_strategy is false.
    // The channel consequently is only populated if skip_reads_length_strategy is false.


    // Initialise channels for outputs as empty to avoid use of conditionals in emit section.
    def ch_seqkit_stats = channel.empty()
    def ch_getseqkitk_kmer = channel.empty()

    if (!params.skip_reads_length_strategy) {
        SEQKIT_STATS(ch_merged_reads)
        GETSEQKITK(SEQKIT_STATS.out.stats)

        // Capture output for emits (avoids using conditionals in emits section)
        ch_seqkit_stats = SEQKIT_STATS.out.stats
        ch_getseqkitk_kmer = GETSEQKITK.out.seqkitkmer_txt

        ch_reads_length_strategy = ch_merged_reads
            .map { meta, reads -> [meta.id, meta, reads] }
            .join(GETSEQKITK.out.seqkitkmer_txt.map { meta, kmer_file ->
                [meta.id, kmer_file.text.trim() as Integer]
            })
            .map { id, meta, reads, seqkit_kmer ->
                // Build k-mer list: add predicted kmer to defaults if valid
                def default_kmers_spades = [21, 33, 55, 77]  // spades recommended defaults
                def max_kmer_spades = 127  // spades max kmer limit

                def default_kmers_megahit = [21, 29, 39, 59, 79, 99, 119, 141] // megahit recommended defaults
                def max_kmer_megahit = 141 // megahit max kmer limit

                // Validate predicted kmer for spades: must be odd, in range [15, 127], and not already in list
                def is_valid_spades = (seqkit_kmer >= 15 &&
                                seqkit_kmer <= max_kmer_spades &&
                                seqkit_kmer % 2 == 1 &&
                                !(seqkit_kmer in default_kmers_spades))

                // Validate predicted kmer for megahit: must be odd, in range [15, 141], and not already in list
                def is_valid_megahit = (seqkit_kmer >= 15 &&
                                seqkit_kmer <= max_kmer_megahit &&
                                seqkit_kmer % 2 == 1 &&
                                !(seqkit_kmer in default_kmers_megahit))

                // Build final k-mer lists for spades
                def kmer_list_spades = is_valid_spades ?
                    (default_kmers_spades + [seqkit_kmer]).sort() :
                    default_kmers_spades

                // Build final k-mer lists for megahit
                def kmer_list_megahit = is_valid_megahit ?
                    (default_kmers_megahit + [seqkit_kmer]).sort() :
                    default_kmers_megahit

                // For assemblers that only take a single k-mer, use the predicted k-mer if it's valid (between 15 and 127), or fall back to a fixed value (25)
                def single_kmer = (seqkit_kmer >= 15 && seqkit_kmer <= 127 && seqkit_kmer % 2 == 1) ? seqkit_kmer : 25

                // Enrich metadata with k-mer strategy and lists
                def enriched_meta = meta + [
                    kmer_strategy: 'reads_length',
                    predicted_kmer: seqkit_kmer,
                    kmer_list_spades: kmer_list_spades.join(','),
                    kmer_list_megahit: kmer_list_megahit.join(','),
                    single_kmer: single_kmer
                ]
                [enriched_meta, reads]
            }
    }
// =================== End of k-mer strategies for genome assembly =======================


// =================== Genome assembly with different assemblers and k-mer strategies =======================

    // Create channel with new meta for downstream processes (after assembly). This channel combines all assemblies from different assemblers and strategies, and maps them to the new meta with updated id.
    def ch_draft_assemblies_input = Channel.empty()

    // ======= Spades assemblies - nested conditionals (assembler × strategy) ======
    // Spades is only run if skip_spades is false. Within that, each strategy is only run if its corresponding skip parameter is false.
    // The channel with Spades assemblies is populated accordingly and mixed into the common ch_draft_assemblies_input channel.
    // Spades needs a tuple with 4 elements as inputs, so we need to map the channel to add empty lists for the other 2 inputs (see PREPROCESSING subworkflow for example)
    // SPADES: [ meta, illumina, pacbio, nanopore ]

    if (!params.skip_spades) {

        if (!params.skip_manual_strategy) {
            // Join manual strategy (merged reads) with unmerged reads
            ch_spades_input_manual = ch_manual_strategy
                .map { meta, merged_reads -> [meta.id, meta, merged_reads] }
                .join(ch_fastp_reads_unmerged.map { meta, pe_reads -> [meta.id, pe_reads] })
                .map { id, meta, merged_reads, pe_reads -> [meta, pe_reads, merged_reads] }

            SPADES_MANUAL(ch_spades_input_manual, [], [])

            // Mix into draft assemblies channel with new meta
            ch_draft_assemblies_input = ch_draft_assemblies_input.mix(
                SPADES_MANUAL.out.scaffolds
                    .map { meta, scaffolds -> createAssemblyMeta(meta, scaffolds, 'spades') }
            )
        }

        if (!params.skip_kmergenie_strategy) {
            // Join kmergenie strategy (merged reads) with unmerged reads
            ch_spades_input_kmergenie = ch_kmergenie_strategy
                .map { meta, merged_reads -> [meta.id, meta, merged_reads] }
                .join(ch_fastp_reads_unmerged.map { meta, pe_reads -> [meta.id, pe_reads] })
                .map { id, meta, merged_reads, pe_reads -> [meta, pe_reads, merged_reads] }

            SPADES_KMERGENIE(ch_spades_input_kmergenie, [], [])

            // Mix into draft assemblies channel with new meta
            ch_draft_assemblies_input = ch_draft_assemblies_input.mix(
                SPADES_KMERGENIE.out.scaffolds
                    .map { meta, scaffolds -> createAssemblyMeta(meta, scaffolds, 'spades') }
            )
        }

        if (!params.skip_reads_length_strategy) {
            // Join reads length strategy (merged reads) with unmerged reads
            ch_spades_input_reads_length = ch_reads_length_strategy
                .map { meta, merged_reads -> [meta.id, meta, merged_reads] }
                .join(ch_fastp_reads_unmerged.map { meta, pe_reads -> [meta.id, pe_reads] })
                .map { id, meta, merged_reads, pe_reads -> [meta, pe_reads, merged_reads] }

            SPADES_READS_LENGTH(ch_spades_input_reads_length, [], [])

            // Mix into draft assemblies channel with new meta
            ch_draft_assemblies_input = ch_draft_assemblies_input.mix(
                SPADES_READS_LENGTH.out.scaffolds
                    .map { meta, scaffolds -> createAssemblyMeta(meta, scaffolds, 'spades') }
            )
        }
    }

    // ======= Megahit assemblies - nested conditionals (assembler × strategy) ======
    // Megahit is only run if skip_megahit is false. Within that, each strategy is only run if its corresponding skip parameter is false.
    // The channel with Megahit assemblies is populated accordingly and mixed into the common ch_draft_assemblies_input channel.
    // Megahit needs a tuple with 3 elements as input.
    // MEGAHIT: [ meta, reads, []]

    if (!params.skip_megahit) {

        if (!params.skip_manual_strategy) {
            ch_megahit_input_manual = ch_manual_strategy.map { meta, reads -> [meta, reads, []] }
            MEGAHIT_MANUAL(ch_megahit_input_manual)

            // Mix into draft assemblies channel with new meta
            ch_draft_assemblies_input = ch_draft_assemblies_input.mix(
                MEGAHIT_MANUAL.out.contigs.map { meta, contigs -> createAssemblyMeta(meta, contigs, 'megahit') }
            )
        }

        if (!params.skip_kmergenie_strategy) {
            ch_megahit_input_kmergenie = ch_kmergenie_strategy.map { meta, reads -> [ meta, reads, [] ] }
            MEGAHIT_KMERGENIE(ch_megahit_input_kmergenie)

            // Mix into draft assemblies channel with new meta
            ch_draft_assemblies_input = ch_draft_assemblies_input.mix(
                MEGAHIT_KMERGENIE.out.contigs.map { meta, contigs -> createAssemblyMeta(meta, contigs, 'megahit') }
            )
        }

        if (!params.skip_reads_length_strategy) {
            ch_megahit_input_reads_length = ch_reads_length_strategy.map { meta, reads -> [ meta, reads, [] ] }
            MEGAHIT_READS_LENGTH(ch_megahit_input_reads_length)

            // Mix into draft assemblies channel with new meta
            ch_draft_assemblies_input = ch_draft_assemblies_input.mix(
                MEGAHIT_READS_LENGTH.out.contigs.map { meta, contigs -> createAssemblyMeta(meta, contigs, 'megahit') }
            )
        }
    }

    // ======= MINIA assemblies - nested conditionals (assembler × strategy) ======
    // Minia is only run if skip_minia is false. Within that, each strategy is only run if its corresponding skip parameter is false.
    // The channel with Minia assemblies is populated accordingly and mixed into the common ch_draft_assemblies_input channel.

    // Initialise channels for outputs as empty to avoid use of conditionals in emit section.
    def ch_minia_contigs_manual = channel.empty()
    def ch_minia_contigs_kmergenie = channel.empty()
    def ch_minia_contigs_reads_length = channel.empty()

    if (!params.skip_minia) {

        if (!params.skip_manual_strategy) {
            MINIA_MANUAL(ch_manual_strategy)

            // Mix into draft assemblies channel with new meta
            ch_draft_assemblies_input = ch_draft_assemblies_input.mix(
                MINIA_MANUAL.out.contigs.map { meta, contigs -> createAssemblyMeta(meta, contigs, 'minia') }
            )
        }

        if (!params.skip_kmergenie_strategy) {
            MINIA_KMERGENIE(ch_kmergenie_strategy)

            // Mix into draft assemblies channel with new meta
            ch_draft_assemblies_input = ch_draft_assemblies_input.mix(
                MINIA_KMERGENIE.out.contigs.map { meta, contigs -> createAssemblyMeta(meta, contigs, 'minia') }
            )
        }

        if (!params.skip_reads_length_strategy) {
            MINIA_READS_LENGTH(ch_reads_length_strategy)

            // Mix into draft assemblies channel with new meta
            ch_draft_assemblies_input = ch_draft_assemblies_input.mix(
                MINIA_READS_LENGTH.out.contigs.map { meta, contigs -> createAssemblyMeta(meta, contigs, 'minia') }
            )
        }
    }

    // ======= ABYSS assemblies - nested conditionals (assembler × strategy) ======
    // ABYSS is only run if skip_abyss is false. Within that, each strategy is only run if its corresponding skip parameter is false.
    // The channel with ABYSS assemblies is populated accordingly and mixed into the common ch_draft_assemblies_input channel.

    if (!params.skip_abyss) {
        if (!params.skip_manual_strategy) {
            // Join manual strategy (merged reads) with original paired-end reads
            ch_abyss_input_manual = ch_manual_strategy
                .map { meta, merged_reads -> [meta.id, meta, merged_reads] }
                .join(ch_fastp_reads.map { meta, pe_reads -> [meta.id, pe_reads] })
                .map { id, meta, merged_reads, pe_reads -> [meta, pe_reads, merged_reads] }
                .multiMap { meta, pe_reads, merged_reads ->
                    input: [ meta, pe_reads, merged_reads ]  // ABYSS expects: [meta, paired_reads, merged_reads]
                    kmer:  params.abyss_kmer
                }

            ABYSS_MANUAL(ch_abyss_input_manual.input, ch_abyss_input_manual.kmer)

            // Mix into draft assemblies channel with new meta
            ch_draft_assemblies_input = ch_draft_assemblies_input.mix(
                ABYSS_MANUAL.out.scaffolds.map { meta, scaffolds -> createAssemblyMeta(meta, scaffolds, 'abyss') }
            )
        }

        if (!params.skip_kmergenie_strategy) {
            // Join kmergenie strategy (merged reads) with original paired-end reads
            ch_abyss_input_kmergenie = ch_kmergenie_strategy
                .map { meta, merged_reads -> [meta.id, meta, merged_reads] }
                .join(ch_fastp_reads.map { meta, pe_reads -> [meta.id, pe_reads] })
                .map { id, meta, merged_reads, pe_reads -> [meta, pe_reads, merged_reads] }
                .multiMap { meta, pe_reads, merged_reads ->
                    input:      [ meta, pe_reads, merged_reads ]  // ABYSS expects: [meta, paired_reads, merged_reads]
                    single_kmer: meta.single_kmer
                }

            ABYSS_KMERGENIE(ch_abyss_input_kmergenie.input, ch_abyss_input_kmergenie.single_kmer)

            // Mix into draft assemblies channel with new meta
            ch_draft_assemblies_input = ch_draft_assemblies_input.mix(
                ABYSS_KMERGENIE.out.scaffolds.map { meta, scaffolds -> createAssemblyMeta(meta, scaffolds, 'abyss') }
            )
        }

        if (!params.skip_reads_length_strategy) {
            // Join reads_length strategy (merged reads) with original paired-end reads
            ch_abyss_input_reads_length = ch_reads_length_strategy
                .map { meta, merged_reads -> [meta.id, meta, merged_reads] }
                .join(ch_fastp_reads.map { meta, pe_reads -> [meta.id, pe_reads] })
                .map { id, meta, merged_reads, pe_reads -> [meta, pe_reads, merged_reads] }
                .multiMap { meta, pe_reads, merged_reads ->
                    input:      [ meta, pe_reads, merged_reads ]  // ABYSS expects: [meta, paired_reads, merged_reads]
                    single_kmer: meta.single_kmer
                }

            ABYSS_READS_LENGTH(ch_abyss_input_reads_length.input, ch_abyss_input_reads_length.single_kmer)

            // Mix into draft assemblies channel with new meta
            ch_draft_assemblies_input = ch_draft_assemblies_input.mix(
                ABYSS_READS_LENGTH.out.scaffolds.map { meta, scaffolds -> createAssemblyMeta(meta, scaffolds, 'abyss') }
            )
        }
    }


    // ======= SPARSEASSEMBLER assemblies - nested conditionals (assembler × strategy) ======
    // SPARSEASSEMBLER is only run if skip_sparseassembler is false. Within that, each strategy is only run if its corresponding skip parameter is false.
    // The channel with SPARSEASSEMBLER assemblies is populated accordingly and mixed into the common ch_draft_assemblies_input channel.
    // SPARSEASSEMBLER is a special case because it can output contigs or scaffolds depending on the parameters used and quality of the reads

    if (!params.skip_sparseassembler) {

        if (!params.skip_manual_strategy) {
            SPARSEASSEMBLER_MANUAL(
                ch_manual_strategy,
                params.sparseassembler_kmer,
                params.sparseassembler_genome_size,
                params.sparseassembler_expected_coverage
            )

            // Mix into draft assemblies channel with new meta
            ch_draft_assemblies_input = ch_draft_assemblies_input.mix(
                SPARSEASSEMBLER_MANUAL.out.scaffolds
                    .concat(SPARSEASSEMBLER_MANUAL.out.contigs)
                    .unique { meta, assembly -> meta.id }
                    .map { meta, assembly -> createAssemblyMeta(meta, assembly, 'sparseassembler') }
            )
        }

        if (!params.skip_kmergenie_strategy) {
            // Create k-mer channel for SPARSEASSEMBLER (needs single kmer value) using multiMap
            ch_sparseassembler_kmergenie = ch_kmergenie_strategy
                .multiMap { meta, reads ->
                    input:       [ meta, reads ]
                    single_kmer: meta.single_kmer
                }

            SPARSEASSEMBLER_KMERGENIE(
                ch_sparseassembler_kmergenie.input,
                ch_sparseassembler_kmergenie.single_kmer,
                params.sparseassembler_genome_size,
                params.sparseassembler_expected_coverage
            )

            // Mix into draft assemblies channel with new meta
            ch_draft_assemblies_input = ch_draft_assemblies_input.mix(
                SPARSEASSEMBLER_KMERGENIE.out.scaffolds
                    .concat(SPARSEASSEMBLER_KMERGENIE.out.contigs)
                    .unique { meta, assembly -> meta.id }
                    .map { meta, assembly -> createAssemblyMeta(meta, assembly, 'sparseassembler') }
            )
        }

        if (!params.skip_reads_length_strategy) {
            // Create k-mer channel for SPARSEASSEMBLER (needs single kmer value) using multiMap
            ch_sparseassembler_reads_length = ch_reads_length_strategy
                .multiMap { meta, reads ->
                    input:       [ meta, reads ]
                    single_kmer: meta.single_kmer
                }

            SPARSEASSEMBLER_READS_LENGTH(
                ch_sparseassembler_reads_length.input,
                ch_sparseassembler_reads_length.single_kmer,
                params.sparseassembler_genome_size,
                params.sparseassembler_expected_coverage
            )

            // Mix into draft assemblies channel with new meta
            ch_draft_assemblies_input = ch_draft_assemblies_input.mix(
                SPARSEASSEMBLER_READS_LENGTH.out.scaffolds
                    .concat(SPARSEASSEMBLER_READS_LENGTH.out.contigs)
                    .unique { meta, assembly -> meta.id }
                    .map { meta, assembly -> createAssemblyMeta(meta, assembly, 'sparseassembler') }
            )
        }
    }

    // ======= Masurca assemblies - nested conditionals (assembler × strategy) ======
    // Masurca is only run if skip_masurca is false. Within that, each strategy is only run if its corresponding skip parameter is false.
    // The channel with Masurca assemblies is populated accordingly and mixed into the common ch_draft_assemblies_input channel.
    // Masurca needs a tuple with 4 elements as inputs, so we need to map the channel to add empty lists for the other 3 inputs (see PREPROCESSING subworkflow for example)
    // MASURCA: [ meta, illumina, jump, pacbio, nanopore ]

    if (!params.skip_masurca) {

        if (!params.skip_manual_strategy) {
            MASURCA_MANUAL(
                ch_manual_strategy.map { meta, reads -> [meta, reads, [], [], []] },
                params.masurca_fragment_mean,
                params.masurca_fragment_stdev,
                0, // no jump reads fragment mean
                0,  // no jump reads fragment stdev
                0, // no extend jump reads
                params.masurca_kmer_size,
                0, // no linking mates
                25, // lhe_coverage - leaving default value for config compatibility, but it's not used as this should be a parameter used for nanopore reads
                0, // mega_reads_one_pass: 0 is default - two passes of mega-reads for slower, but higher quality assembly
                300, // limit_jump_coverage - leaving default value for config compatibility, but it's not used as this should be a parameter used for jump reads
                params.masurca_ca_parameters, // cgwErrorRate
                1, // do attempt to close gaps (we can leave this hardcoded)
                params.masurca_jf_size
            )

            // Mix into draft assemblies channel with new meta
            ch_draft_assemblies_input = ch_draft_assemblies_input.mix(
                MASURCA_MANUAL.out.scaffolds
                    .map { meta, scaffolds -> createAssemblyMeta(meta, scaffolds, 'masurca') }
            )
        }

        if (!params.skip_kmergenie_strategy) {

            // Create k-mer channel for Masurca (needs single kmer value) using multiMap
            ch_masurca_kmergenie = ch_kmergenie_strategy
                .multiMap { meta, reads ->
                    input:      [ meta, reads, [], [], [] ]
                    single_kmer: meta.single_kmer
                }

            MASURCA_KMERGENIE(
                ch_masurca_kmergenie.input,
                params.masurca_fragment_mean,
                params.masurca_fragment_stdev,
                0, // no jump reads fragment mean
                0,  // no jump reads fragment stdev
                0, // no extend jump reads
                ch_masurca_kmergenie.single_kmer,
                0, // no linking mates
                25, // lhe_coverage - leaving default value for config compatibility, but it's not used as this should be a parameter used for nanopore reads
                0, // mega_reads_one_pass: 0 is default - two passes of mega-reads for slower, but higher quality assembly
                300, // limit_jump_coverage - leaving default value for config compatibility, but it's not used as this should be a parameter used for jump reads
                params.masurca_ca_parameters, // cgwErrorRate
                1, // do attempt to close gaps (we can leave this hardcoded)
                params.masurca_jf_size
            )

            // Mix into draft assemblies channel with new meta
            ch_draft_assemblies_input = ch_draft_assemblies_input.mix(
                MASURCA_KMERGENIE.out.scaffolds
                    .map { meta, scaffolds -> createAssemblyMeta(meta, scaffolds, 'masurca') }
            )
        }

        if (!params.skip_reads_length_strategy) {

            // Create k-mer channel for Masurca (needs single kmer value) using multiMap
            ch_masurca_reads_length = ch_reads_length_strategy
                .multiMap { meta, reads ->
                    input:      [ meta, reads, [], [], [] ]
                    single_kmer: meta.single_kmer
                }

            MASURCA_READS_LENGTH(
                ch_masurca_reads_length.input,
                params.masurca_fragment_mean,
                params.masurca_fragment_stdev,
                0, // no jump reads fragment mean
                0,  // no jump reads fragment stdev
                0, // no extend jump reads
                ch_masurca_reads_length.single_kmer,
                0, // no linking mates
                25, // lhe_coverage - leaving default value for config compatibility, but it's not used as this should be a parameter used for nanopore reads
                0, // mega_reads_one_pass: 0 is default - two passes of mega-reads for slower, but higher quality assembly
                300, // limit_jump_coverage - leaving default value for config compatibility, but it's not used as this should be a parameter used for jump reads
                params.masurca_ca_parameters, // cgwErrorRate
                1, // do attempt to close gaps (we can leave this hardcoded)
                params.masurca_jf_size
            )

            // Mix into draft assemblies channel with new meta
            ch_draft_assemblies_input = ch_draft_assemblies_input.mix(
                MASURCA_READS_LENGTH.out.scaffolds
                    .map { meta, scaffolds -> createAssemblyMeta(meta, scaffolds, 'masurca') }
            )
        }
    }

    emit:
    // K-mer strategy outputs - might be useful to collect results downstream
    seqkit_stats                                = ch_seqkit_stats
    getseqkitk_kmer                             = ch_getseqkitk_kmer
    getkmergeniek_k                             = ch_getkmergeniek_k
    // Assemblies - all assemblers and strategies mixed into one channel
    draft_assemblies_merged = ch_draft_assemblies_input
}
