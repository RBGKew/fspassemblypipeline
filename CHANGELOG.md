# nf-core/fspassemblypipeline: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v1.0.0dev - [unreleased<!-- TODO nf-core: replace with date on release -->]

Initial release of nf-core/fspassemblypipeline, created with the [nf-core](https://nf-co.re/) template.

### `Added`

- 09/02/2026 - Initialised GENOME_ASSEMBLY subworkflow.
- 09/02/2026 - Added SEQKIT_STATS module.
- 10/02/2026 - Added FASTK_FASTK module.
- 10/02/2026 - Added SPADES module.
- 10/02/2026 - Added MEGAHIT module.
- 10/02/2026 - Added MINIA module.
- 12/02/2026 - Added RENAME_ASSEMBLIES local module.
- 12/02/2026 - Added BUSCO_BUSCO module.
- 12/02/2026 - Added MERQURYFK_MERQURYFK module.
- 12/02/2026 - Added QUAST module.
- 25/02/2026 - Added kmergenie local module.
- 25/02/2026 - Added getkmergeniek local module.
- 27/02/2026 - Created nf-core module for kmergenie.
- 03/03/2026 - Added kmergenie to the pipeline as nf-core module.
- 10/03/2026 - Added GENOMESCOPE2 module.
- 10/03/2026 - Added GENESCOPEFK module.
- 10/03/2026 - Added KMC module.
- 10/03/2026 - Added FALCO_QCSTAT_COMPILING local module.
- 10/03/2026 - Added FQSTAT local module.
- 10/03/2026 - Added FQSTATSUMMARY local module.
- 10/03/2026 - Added KMER_STAT_SUMMARY local module.
- 10/03/2026 - Added PREPROCESSING subworkflow skeleton (main.nf and meta.yml) and supporting scripts in bin.
- 11/03/2026 - Added FASTK_HISTEX module.
- 17/03/2026 - Added sparseassembler as local module.
- 17/03/2026 - Added ABYSS_ABYSSPE.
- 17/03/2026 - Added FASTP parameter settings in nextflow.config to nf-core-style individual parameters.
- 17/03/2026 - Added FASTP module argument construction to conf/modules.config and added FASTP-specific publishDir output path.
- 26/03/2026 - Added nf-core-style configurable FastK/Histex parameters (`-v`, `-p`, `fastk_extra_args`, `-G`, `histex_extra_args`) to nextflow.config and nextflow_schema.json.
- 26/03/2026 - Added dedicated `FALCO_AFTER_MERGE` preprocessing invocation and per-sample Falco publish layout under `falco/<sample>/{raw,trimmed,merge}`.
- 27/03/2026 - Brought back missing `MULTIQC` module include in `workflows/fspassemblypipeline.nf` as it causes error I don't want to fix.
- 02/04/2026 - Synced changelog entries from `Wu_fspassemblypipeline` and documented 02/04 debugging updates in this repository changelog.
- 16/04/2026 - Implemented selection of busco lineage per sample
- 21/04/2026 - implemented kmergenie and reads_length kmer strategies
- 29/04/2026 - Kmer strategies (manual, kmergenie, reads_length) implemented. Now the assemblers are run x3 with different kmer sizes.
- 08/05/2026 - The workflow can now be run selecting some assemblers and kmer strategies, no need to run them all if it's not necessary.
- 08/05/2026 - If no taxonomy is provided, only the preprocessing will be run.
- 22/05/2026 - Thw workflow now runs using pre-processed paired-end reads (GENOME_ASSEMBLY), and merged reads (GENOME_ASSEMBLY_MERGED). It possible select one or the other with the parameters `use_paired_reads` and `use_merged_reads`.
- 11/06/2026 - The QC of the draft assemblies has been moved to a separate subworkflow (SELECT_BEST_ASSEMBLY_AND_QC). This subworkflow collects the assemblies from GENOME_ASSEMBLY and GENOME_ASSEMBLY_MERGED, runs BUSCO, QUAST and Merqury for all of them, and then runs a script that selects the best assembly based on single-copy complete BUSCOs obtained using the specific lineage. If two or more assemblies have the same BUSCO score, the auN calculated by QUAST is used as a metric to evaluate the contiguity, and the most contiguous assembly is selected.
- 17/06/2026 - SELECT_BEST_ASSEMBLY_AND_QC now runs pypolca to improve the selected best assembly (error correction). The corrected assembly is then QCed using BUSCO, QUAST and Merqury, and bwa-mem2 and samtools are used to analyse the coverage. A coverage visualisation script is also run for rapid analysis.
- 18/06/2026 - Added Masurca assembler.
- 23/06/2026 - Added local module for spades using merged reads.

### `Fixed`

- 12/02/2026 - Renamed assemblies to avoid conflicts in downstream modules.
- 12/02/2026 - New output directory structure.
- 12/02/2026 - User can set extra params from `nextflow.config`.
- 10/03/2026 - Removed obsolete script `bin/draw_distinct_kmer_vs_reads.py` for preprocessing.
- 11/03/2026 - Updated preprocessing and FQSTATSUMMARY inputs using def function to pass file paths explicitly for nf-core compliance.
- 11/03/2026 - Updated FALCO_QCSTAT_COMPILING parsing to process Falco outputs by file type and simplified environment dependencies.
- 11/03/2026 - Removed KMC module from the pipeline.
- 12/03/2026 - Added FASTK log output and adjusted GENESCOPEFK output folder/suffix for downstream compatibility.
- 13/03/2026 - Fixed kmergenie nf-core module (missing log).
- 13/03/2026 - Updated FASTK and GENESCOPEFK nf-core modules to emit log files and updated the modules from nf-core rather than making local changes.
- 16/03/2026 - Updated GENOME_ASSEMBLY subworkflow to incorporate both GENESCOPEFK peak 1 and peak 2.
- 17/03/2026 - Removed tests and meta.yml from getkmergeniek local module as it's not needed.
- 17/03/2026 - Aligned FASTK_HISTEX outputs to `*.hist.txt` and updated downstream preprocessing outputs/metadata and k-mer summary histogram parsing accordingly.
- 26/03/2026 - Updated FASTK_FASTK/FASTK_HISTEX module argument construction in conf/modules.config and set FastK `-p`/`-v` to be enabled by default to preserve required outputs.
- 26/03/2026 - Updated falco to new topic version.
- 26/03/2026 - Updated Falco QC statistics aggregation to include `raw`, `trimmed`, and `merge`, and renamed outputs to `QC_raw_result.txt`, `QC_trimmed_result.txt`, and `QC_merge_result.txt`.
- 26/03/2026 - Updated `FQSTAT` to analyse trimmed R1/R2 plus merged reads; fixed `FQSTAT_SUMMARY` parsing of `Total:`/`Average:` fields.
- 27/03/2026 - Updated `GENESCOPEFK_P1`/`GENESCOPEFK_P2` publish paths and runtime args to use separated `p1`/`p2` outputs with `-p 1 -k 17` and `-p 2 -k 17`.
- 27/03/2026 - Simplified `KMER_STAT_SUMMARY` input staging and set `statistics_all.csv` `peak_positions` from `peak_1` `genescopefk.log` `kcov` while keeping the original three Python script execution order.
- 27/03/2026 - Synced local `fastk/histex` and `genescopefk` nf-core module copies and aligned schema defaults for `fastp_extra_args`, `fastk_extra_args`, and `histex_extra_args` with lint expectations.
- 31/03/2026 - Moved the large raw reads and bam to large file system (lfs), and stored as lfs pointers.
- 01/04/2026 - Split the preprocessing PR to four parts, each part origins from newly forked dev from kew repository.
- 02/04/2026 - Merged local modules and samplesheet PR, updating the part of local modules.
- 10/04/2026 - updated samplesheet and updated the following modules to latest version and use topic channels: busco, fastp, megahit, minia, quast, seqkit, spades, multiqc.

### `Dependencies`

### `Deprecated`

- 13/03/2026 - Discarded the changes from "fixed" 12/03/2026.
