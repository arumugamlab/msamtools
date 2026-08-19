# msamtools: microbiome-related extension to samtools

**msamtools** provides functions commonly used when processing SAM/BAM
alignments in microbiome data analysis, particularly shotgun metagenomics and
metatranscriptomics.

It is designed to work with alignment streams in a similar spirit to
**samtools**, allowing filtering, profiling, coverage estimation, and alignment
summarization to be combined in pipelines without creating unnecessary
intermediate files.

## Commands

msamtools provides four main commands:

- [`filter`](filter.md) — filter alignments based on alignment statistics
- [`profile`](profile.md) — estimate abundance profiles of reference sequences or genomes
- [`coverage`](coverage.md) — estimate per-position or per-sequence coverage
- [`summary`](summary.md) — summarize alignment statistics

See the [example workflows](examples.md) for ways to combine msamtools commands
with each other and with samtools.

## Getting started

Start with:

- [Installation](installation.md)
- [Input requirements](input-requirements.md)
- [Compatibility](compatibility.md)

For reproducible validation of the profiling functionality, see the
[profile validation report](validation/profile_validation.md).
