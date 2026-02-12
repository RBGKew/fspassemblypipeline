# nf-core/fsptest: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v1.0.0dev - [date]

Initial release of nf-core/fsptest, created with the [nf-core](https://nf-co.re/) template.

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

### `Fixed`

09/02/2026 - Commented out ch_versions in fsptest.nf, as it was causing the failure of the test.
12/02/2026 - Renamed assemblies to avoid conflicts in downstream modules.
12/02/2026 - New output directory structure.
12/02/2026 - User cans set extra params from `nextflow.config`.

### `Dependencies`

### `Deprecated`
