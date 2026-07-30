include { CREATE_PROJECT_YAML           } from '../../../modules/local/createyml/main'
include { BLOBTOOLKIT_CREATEBLOBDIR     } from '../../../modules/local/blobtoolkit_create/main'
include { SAMTOOLS_INDEX as SAMTOOLS_CSI} from '../../../modules/nf-core/samtools/index/main'


workflow BLOBTOOLS {

    take:
    ch_samplesheet
    ch_blobtools_taxonomy

    main:
    // a channel for each input is create and finally all channels combined into one tuple to give blobtools a single tuple.
    // This is to allow for easier adding subtracting of inputs for blobtools but just removing or adding channels. Also blobtools wants a single tuple as input


        // Create YAML files from samplesheet metadata
        ch_yaml_input = ch_samplesheet
            .map { meta, files ->
                tuple(meta, files[0])  // meta and fasta file
            }

        CREATE_PROJECT_YAML(ch_yaml_input)

        // Debug: Check what YAML output looks like
        CREATE_PROJECT_YAML.out.yaml.view { meta, yaml ->
            "YAML created: ${meta.id} -> ${yaml}"
        }

        // Prepare channels for joining
        ch_samplesheet_keyed = ch_samplesheet
            .map { meta, files ->
                tuple(meta.id, meta, files[0], files[1], files[2])  // fasta and bam
            }
            .view { "Samplesheet keyed: ${it}" }

        ch_yaml_keyed = CREATE_PROJECT_YAML.out.yaml
            .map { meta, yaml ->
                tuple(meta.id, yaml)
            }
            .view { "YAML keyed: ${it}" }

        ch_taxonomy_keyed = ch_blobtools_taxonomy
            .map { meta, file ->
                tuple(meta.id, file)
            }


        ch_samtools = ch_samplesheet
            .map { meta, files ->
                tuple(meta, files[1])  // fasta and bam
            }

        SAMTOOLS_CSI(ch_samtools)


        ch_samtools_keyed = SAMTOOLS_CSI.out.index
            .map { meta, index ->
                tuple(meta.id, index)
            }
            .view { "index keyed: ${it}" }


        // Join and combine with BUSCO
        ch_btk = ch_samplesheet_keyed
            .join(ch_yaml_keyed)
            .join(ch_samtools_keyed)
            .join(ch_taxonomy_keyed)
            .map { id, meta, fasta, bam, busco, yaml, index, taxonomy ->
                tuple(meta, fasta, bam, busco, yaml, index, taxonomy)
            }
            .view { "Final input to BLOBTOOLKIT (meta, fasta, bam, busco, yaml, index, taxonomy): ${it}" }
        //
        // Create Blobtools dataset files
        //
        BLOBTOOLKIT_CREATEBLOBDIR(ch_btk)

    emit:
    blobdir  = BLOBTOOLKIT_CREATEBLOBDIR.out.blobdir
}
