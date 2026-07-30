include { TIARA_TIARA        } from '../../../modules/nf-core/tiara/tiara/main'
include { FCSGX_RUNGX        } from '../../../modules/nf-core/fcsgx/rungx/main'
include { CONVERTFCSRPT      } from '../../../modules/local/convertfcsrpt/main'
include { COMPARISON         } from '../../../modules/local/comparison/main'

workflow CONTAMINATION_DETECTION {

    take:
    ch_samplesheet

    main:
    ch_tiara = ch_samplesheet
        .map { meta, files ->
            tuple(
                [id: meta.id],      // Create new meta with just id
                files[0]            // First file in the list is the fasta
            )
        }

    ch_fcsgx = ch_samplesheet
        .map { meta, files ->
            tuple(
                [id: meta.id],      // Create new meta with just id
                meta.taxid,         // Extract taxid from meta
                files[0]            // First file in the list is the fasta
            )
        }

    // Run Tiara classification
    TIARA_TIARA(ch_tiara)

    // Run FCS-GX contamination screening
    FCSGX_RUNGX(ch_fcsgx, params.db_path, params.ramdisk_path ?:[])

    // Convert FCS-GX report format
    CONVERTFCSRPT(FCSGX_RUNGX.out.taxonomy_report)

    // Compare Tiara and FCS-GX results
    ch_comparison_input = TIARA_TIARA.out.classifications
        .join(CONVERTFCSRPT.out.fcs_report_reformatted)

    COMPARISON(ch_comparison_input)

    emit:
    tiara_classifications  = TIARA_TIARA.out.classifications
    taxonomy_report        = FCSGX_RUNGX.out.taxonomy_report
    fcsgx_reformatted      = CONVERTFCSRPT.out.fcs_report_reformatted
    blobtools_taxonomy     = COMPARISON.out.blobtools_taxonomy
}
