# nf-core/fspassemblypipeline: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v1.0.0dev - [date]

Initial release of nf-core/fspassemblypipeline, created with the [nf-core](https://nf-co.re/) template.

### `Added`

09/02/2026 - Initialised GENOME_ASSEMBLY subworkflow.
09/02/2026 - Added SEQKIT_STATS module.
10/02/2026 - Added FASTK_FASTK module.
10/02/2026 - Added SPADES module.
10/02/2026 - Added MEGAHIT module.
10/02/2026 - Added MINIA module.
12/02/2026 - Added RENAME_ASSEMBLIES local module.
12/02/2026 - Added BUSCO_BUSCO module.
12/02/2026 - Added MERQURYFK_MERQURYFK module.
12/02/2026 - Added QUAST module.
10/03/2026 - Added GENOMESCOPE2 module.
10/03/2026 - Added GENESCOPEFK module.
10/03/2026 - Added KMC module.
10/03/2026 - Added FALCO_QCSTAT_COMPILING local module.
10/03/2026 - Added FQSTAT local module.
10/03/2026 - Added FQSTATSUMMARY local module.
10/03/2026 - Added KMER_STAT_SUMMARY local module.
10/03/2026 - Added PREPROCESSING subworkflow skeleton (main.nf and meta.yml) and supporting scripts in bin.
11/03/2026 - Added FASTK_HISTEX module.
17/03/2026 - Added FASTP parameter settings in nextflow.config to nf-core-style individual parameters.
17/03/2026 - Added FASTP module argument construction to conf/modules.config and added FASTP-specific publishDir output path.
26/03/2026 - Added nf-core-style configurable FastK/Histex parameters (`-v`, `-p`, `fastk_extra_args`, `-G`, `histex_extra_args`) to nextflow.config and nextflow_schema.json.
26/03/2026 - Added dedicated `FALCO_AFTER_MERGE` preprocessing invocation and per-sample Falco publish layout under `falco/<sample>/{raw,trimmed,merge}`.
27/03/2026 - Brought back missing `MULTIQC` module include in `workflows/fspassemblypipeline.nf` as it causes error I don't want to fix.

### `Fixed`

09/02/2026 - Commented out ch_versions in fsptest.nf, as it was causing the failure of the test.
12/02/2026 - Renamed assemblies to avoid conflicts in downstream modules.
12/02/2026 - New output directory structure.
12/02/2026 - User cans set extra params from `nextflow.config`.
10/03/2026 - Updated KMC module version/input wiring and fixed lint-related configuration.
10/03/2026 - Removed obsolete script bin/draw_distinct_kmer_vs_reads.py for preprocessing.
11/03/2026 - Updated preprocessing and FQSTATSUMMARY inputs using def function to pass file paths explicitly for nf-core compliance.
11/03/2026 - Updated FALCO_QCSTAT_COMPILING parsing to process Falco outputs by file type and simplified environment dependencies.
11/03/2026 - Removed KMC module from the pipeline.
12/03/2026 - Added FASTK log output and adjusted GENESCOPEFK output folder/suffix for downstream compatibility.
13/03/2026 - Updated FASTK and GENESCOPEFK nf-core modules to emit log files and updated the modules from nf-core rather than making local changes.
16/03/2026 - Updated GENOME_ASSEMBLY subworkflow to incorporate both GENESCOPEFK peak 1 and peak 2.
17/03/2026 - Aligned FASTK_HISTEX outputs to `*.hist.txt` and updated downstream preprocessing outputs/metadata and k-mer summary histogram parsing accordingly.
26/03/2026 - Updated FASTK_FASTK/FASTK_HISTEX module argument construction in conf/modules.config and set FastK `-p`/`-v` to be enabled by default to preserve required outputs.
26/03/2026 - Updated falco to new topic version.
26/03/2026 - Updated Falco QC statistics aggregation to include `raw`, `trimmed`, and `merge`, and renamed outputs to `QC_raw_result.txt`, `QC_trimmed_result.txt`, and `QC_merge_result.txt`.
26/03/2026 - Updated `FQSTAT` to analyse trimmed R1/R2 plus merged reads; fixed `FQSTAT_SUMMARY` parsing of `Total:`/`Average:` fields.
27/03/2026 - Updated `GENESCOPEFK_P1`/`GENESCOPEFK_P2` publish paths and runtime args to use separated `p1`/`p2` outputs with `-p 1 -k 17` and `-p 2 -k 17`.
27/03/2026 - Simplified `KMER_STAT_SUMMARY` input staging and set `statistics_all.csv` `peak_positions` from `peak_1` `genescopefk.log` `kcov` while keeping the original three Python script execution order.
27/03/2026 - Synced local `fastk/histex` and `genescopefk` nf-core module copies and aligned schema defaults for `fastp_extra_args`, `fastk_extra_args`, and `histex_extra_args` with lint expectations.

### `Dependencies`

### `Deprecated`

13/03/2026 - Discarded the changes from "fixed" 12/03/2026.
