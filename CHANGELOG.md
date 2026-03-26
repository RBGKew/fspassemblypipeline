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

### `Dependencies`

### `Deprecated`

13/03/2026 - Discarded the changes from "fixed" 12/03/2026.
