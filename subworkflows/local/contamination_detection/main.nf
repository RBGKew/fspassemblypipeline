include { TIARA_TIARA        } from '../../../modules/nf-core/tiara/tiara/main'
include { FCSGX_RUNGX        } from '../../../modules/nf-core/fcsgx/rungx/main' 
include { CONVERTFCSRPT      } from '../../../modules/local/convertfcsrpt/main'
include { COMPARISON         } from '../../../modules/local/comparison/main'

workflow CONTAMINATION_DETECTION {

    take:
    ch_assemblies       // channel: [val(meta), path(assembly)] - genome assemblies (FASTA or BAM files)
    ch_ramdisk_path     // value: ramdisk path or empty list
    ch_db_path          // value: database path
   
    main:
    // Run Tiara classification
    TIARA_TIARA(ch_assemblies)

    // Prepare input for FCS-GX: use per-sample taxon_id from metadata, fall back to global params.taxid
    ch_fcs_gx = ch_assemblies.map { meta, assembly -> 
        def taxid = meta.taxon_id ?: params.taxid
        [meta, taxid, assembly]
    }

    // Run FCS-GX contamination screening
    FCSGX_RUNGX(ch_fcs_gx, ch_db_path, ch_ramdisk_path)

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
    
    versions = TIARA_TIARA.out.versions.mix(FCSGX_RUNGX.out.versions) 
}



