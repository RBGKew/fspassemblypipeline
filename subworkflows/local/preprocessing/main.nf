// Subworkflow implementation for read preprocessing and QC/statistics collation.

include { FALCO as FALCO_RAW            } from '../../../modules/nf-core/falco/main'
include { FALCO as FALCO_AFTER_FASTP    } from '../../../modules/nf-core/falco/main'
include { FALCO as FALCO_AFTER_MERGE    } from '../../../modules/nf-core/falco/main'
include { FASTP as FASTP_TRIM           } from '../../../modules/nf-core/fastp/main'
include { FASTP as FASTP_MERGE          } from '../../../modules/nf-core/fastp/main'
include { FASTK_FASTK                   } from '../../../modules/nf-core/fastk/fastk/main'
include { FASTK_HISTEX                  } from '../../../modules/nf-core/fastk/histex/main'
include { GENESCOPEFK as GENESCOPEFK_P1 } from '../../../modules/nf-core/genescopefk/main'
include { GENESCOPEFK as GENESCOPEFK_P2 } from '../../../modules/nf-core/genescopefk/main'
include { FALCO_QCSTAT_COMPILING        } from '../../../modules/local/falco_qcstat_compiling/main'
include { FQSTAT                        } from '../../../modules/local/fqstat/main'
include { FQSTAT_SUMMARY                } from '../../../modules/local/fqstat_summary/main'
include { KMER_STAT_SUMMARY             } from '../../../modules/local/kmer_stat_summary/main'

workflow PREPROCESSING {

    take:
    ch_samplesheet // channel: [ val(meta), [ reads ] ]

    main:
    ch_versions = Channel.empty()

    FALCO_RAW (
        ch_samplesheet
    )
    ch_versions = ch_versions.mix( FALCO_RAW.out.versions_falco )

    // FASTP expects [meta, reads, adapter_fasta].
    def fastp_adapter_fasta = params.fastp_adapter_fasta ? file(params.fastp_adapter_fasta, checkIfExists: true) : []

    // Step 1: quality trimming and filtering, keeping full trimmed R1/R2 output.
    ch_samplesheet_fastp_trim = ch_samplesheet.map { meta, reads -> [ meta, reads, fastp_adapter_fasta ] }

    FASTP_TRIM (
        ch_samplesheet_fastp_trim,
        false,
        false,
        false
    )
    ch_versions = ch_versions.mix( FASTP_TRIM.out.versions_fastp )

    // Step 2: merge from trimmed reads while retaining unmerged and merged outputs.
    ch_samplesheet_fastp_merge = FASTP_TRIM.out.reads.map { meta, reads -> [ meta, reads, [] ] }

    FASTP_MERGE (
        ch_samplesheet_fastp_merge,
        false,
        false,
        true
    )
    ch_versions = ch_versions.mix( FASTP_MERGE.out.versions_fastp )

    FALCO_AFTER_FASTP (
        FASTP_TRIM.out.reads
    )
    ch_versions = ch_versions.mix( FALCO_AFTER_FASTP.out.versions_falco )

    FALCO_AFTER_MERGE (
        FASTP_MERGE.out.reads_merged
    )
    ch_versions = ch_versions.mix( FALCO_AFTER_MERGE.out.versions_falco )

    ch_falco_qc_compiling_input = FALCO_RAW.out.txt
        .flatMap { meta, txt_files ->
            def files = txt_files instanceof List ? txt_files : [txt_files]
            files.collect { txt -> [ 'raw', txt ] }
        }
        .mix(
            FALCO_AFTER_FASTP.out.txt
                .flatMap { meta, txt_files ->
                    def files = txt_files instanceof List ? txt_files : [txt_files]
                    files.collect { txt -> [ 'trimmed', txt ] }
                }
        )
        .mix(
            FALCO_AFTER_MERGE.out.txt
                .flatMap { meta, txt_files ->
                    def files = txt_files instanceof List ? txt_files : [txt_files]
                    files.collect { txt -> [ 'merge', txt ] }
                }
        )
        .groupTuple(by: 0)

    FALCO_QCSTAT_COMPILING (
        ch_falco_qc_compiling_input
    )

// make a mixed channel for trimmed R1/R2 and merged reads to be processed by FQSTAT.
    ch_fqstat_input = FASTP_TRIM.out.reads
        .mix( FASTP_MERGE.out.reads_merged )

    FQSTAT (
        ch_fqstat_input
    )

    ch_fqstat_summary_input = FQSTAT.out.stats
        .map { meta, stats_file -> stats_file }
        .flatten()
        .collect()
        .filter { it }

    FQSTAT_SUMMARY (
        ch_fqstat_summary_input
    )
    FASTK_FASTK (
        FASTP_TRIM.out.reads
    )
    ch_versions = ch_versions.mix( FASTK_FASTK.out.versions_fastk )

    FASTK_HISTEX (
        FASTK_FASTK.out.hist
    )
    ch_versions = ch_versions.mix( FASTK_HISTEX.out.versions_fastk )

    GENESCOPEFK_P1 (
        FASTK_HISTEX.out.hist
    )
    ch_versions = ch_versions.mix( GENESCOPEFK_P1.out.versions_genescopefk )

    GENESCOPEFK_P2 (
        FASTK_HISTEX.out.hist
    )
    ch_versions = ch_versions.mix( GENESCOPEFK_P2.out.versions_genescopefk )
// Prepare a mixed channel of FastK histograms and Genomescope summaries for KMER_STAT_SUMMARY.
    ch_fastk_hist_input = FASTK_HISTEX.out.hist
        .map { meta, hist -> hist }
        .collect()
        .filter { it }

    ch_fastk_log_input = FASTK_FASTK.out.log
        .map { meta, log -> log }
        .collect()
        .filter { it }

    ch_genomescope_p1_summary_input = GENESCOPEFK_P1.out.summary
        .map { meta, summary -> summary }
        .collect()
        .filter { it }

    ch_genomescope_p2_summary_input = GENESCOPEFK_P2.out.summary
        .map { meta, summary -> summary }
        .collect()
        .filter { it }

    ch_genomescope_p1_log_input = GENESCOPEFK_P1.out.log
        .map { meta, log -> log }
        .collect()
        .filter { it }

    ch_genomescope_p2_log_input = GENESCOPEFK_P2.out.log
        .map { meta, log -> log }
        .collect()
        .filter { it }

    KMER_STAT_SUMMARY (
        ch_fastk_hist_input,
        ch_fastk_log_input,
        ch_genomescope_p1_summary_input,
        ch_genomescope_p2_summary_input,
        ch_genomescope_p1_log_input,
        ch_genomescope_p2_log_input
    )

    emit:
    fastp_reads            = FASTP_TRIM.out.reads // complete trimmed R1 && R2
    fastp_reads_merged     = FASTP_MERGE.out.reads_merged
    falco_raw_html         = FALCO_RAW.out.html
    falco_after_fastp_html = FALCO_AFTER_FASTP.out.html
    falco_after_merge_html = FALCO_AFTER_MERGE.out.html
    falco_qc_stats         = FALCO_QCSTAT_COMPILING.out.stats
    fq_stats               = FQSTAT.out.stats
    fq_stats_summary       = FQSTAT_SUMMARY.out.summary
    fastk_ktab             = FASTK_FASTK.out.ktab
    fastk_hist             = FASTK_FASTK.out.hist
    histex_txt             = FASTK_HISTEX.out.hist
    genomescope_summary    = GENESCOPEFK_P1.out.summary.mix( GENESCOPEFK_P2.out.summary )
    kmer_stats_summary     = KMER_STAT_SUMMARY.out.summary
    versions               = ch_versions
}
