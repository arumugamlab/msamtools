# msamtools: microbiome-related extension to samtools

**msamtools** provides functions commonly used when processing SAM/BAM alignments in microbiome data analysis, particularly shotgun metagenomics and metatranscriptomics.

It is designed to work with alignment streams in a similar spirit to **samtools**, allowing filtering, profiling, coverage estimation, and alignment summarization to be combined in pipelines without creating unnecessary intermediate files.

**Documentation:** [https://msamtools.readthedocs.io/](https://msamtools.readthedocs.io/)

## Installation

The easiest way to install **msamtools** is through [Bioconda](https://bioconda.github.io/):

```bash
conda install -c bioconda msamtools
```

To create a separate environment:

```bash
conda create -n msamtools -c conda-forge -c bioconda msamtools
conda activate msamtools
```

Source builds and container-based usage are described in the [installation documentation](docs/installation.md).

## Commands

**msamtools** provides four commands:

- **`filter`** — filter alignments based on alignment statistics
- **`profile`** — estimate abundance profiles of reference sequences or genomes
- **`coverage`** — estimate per-base or per-sequence read coverage
- **`summary`** — summarize alignment statistics

Run:

```bash
msamtools help
```

for the command overview, or:

```bash
msamtools <command> --help
```

for command-specific options.

## Example

Commands can be streamed together. For example:

```bash
msamtools filter -b -u -l 80 -p 95 -z 80 --besthit sample.bam \
    | msamtools profile --multi=proportional --label=sample --unit=rel -o sample.profile.txt.gz -
```

More workflows are available in the [examples documentation](docs/examples.md).

## Documentation

Detailed documentation is available under [`docs/`](docs/index.md):

- [Installation](docs/installation.md)
- [Examples](docs/examples.md)
- [`filter`](docs/filter.md)
- [`profile`](docs/profile.md)
- [`coverage`](docs/coverage.md)
- [`summary`](docs/summary.md)
- [Profile validation](docs/validation/profile_validation.md)

The documentation is built using [MkDocs](https://www.mkdocs.org/).

To preview it locally:

```bash
mkdocs serve
```

## Citation

If you use **msamtools** in your work, please see [`CITATION.cff`](CITATION.cff) for citation information.

## License

**msamtools** is distributed under the terms described in [`LICENSE`](LICENSE).
