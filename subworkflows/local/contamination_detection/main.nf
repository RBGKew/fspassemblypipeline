// TODO nf-core: If in doubt look at other nf-core/subworkflows to see how we are doing things! :)
//               https://github.com/nf-core/modules/tree/master/subworkflows
//               You can also ask for help via your pull request or on the #subworkflows channel on the nf-core Slack workspace:
//               https://nf-co.re/join
// TODO nf-core: A subworkflow SHOULD import at least two modules

//include { SAMTOOLS_SORT      } from '../../../modules/nf-core/samtools/sort/main'
//include { SAMTOOLS_INDEX     } from '../../../modules/nf-core/samtools/index/main'
include   { TIARA_TIARA        } from '../../../modules/nf-core/tiara/tiara/main'

workflow CONTAMINATION_DETECTION {

    take:
    // TODO nf-core: edit input (take) channels
    ch_tiara_input // [val(meta), path(assemblies)]

    main:
    TIARA_TIARA(ch_tiara_input)  // call module

    emit:
    classifications  = TIARA_TIARA.out.classifications  
    versions = TIARA_TIARA.out.versions 
    
    
}




