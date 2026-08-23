# Installation

## System requirements

**msamtools** is intended for Linux and other UNIX-like systems. A macOS build
is also available through Bioconda.

**msamtools** uses HTSlib for reading and writing SAM/BAM files. Source builds
download and build the HTSlib version specified by the release. The current
development version uses HTSlib 1.24.

## Recommended installation using Conda

The easiest way to install **msamtools** and its required dependencies is
through the Bioconda channel.

If you are already within a Conda environment:

```bash
conda install -c bioconda msamtools
msamtools help
```

To create a separate environment:

```bash
conda create -n msamtools -c conda-forge -c bioconda msamtools
conda activate msamtools
msamtools help
```

If you also need the **samtools** executable for your analysis, install it in
the same environment:

```bash
conda create -n msamtools -c conda-forge -c bioconda msamtools samtools
conda activate msamtools
msamtools help
samtools help
```

If you cannot install through Conda, see
[Advanced installation](#advanced-installation).

## Using container images without installing locally

**msamtools** is available in container images that can be used, for example,
in workflow systems such as Snakemake.

Two container options are available.

### Bioconda container

The Bioconda release provides a corresponding container image containing
**msamtools**.

For example:

```snakemake
singularity: 'docker://quay.io/biocontainers/msamtools:MSAM_VERSION--BIOCONDA_BUILD'
```

The image can then be used by Snakemake with container support enabled.

### msamtools + samtools container

For workflows that pipe between **msamtools** and **samtools**, we also provide
a container containing both programs.

The samtools version in this image is built to match the HTSlib version used by
msamtools. For the current development version, this is samtools/HTSlib 1.24.

For example, if a BAM file is coordinate-sorted and needs to be name-sorted
before profiling:

```snakemake
rule profile_sample:
    input: "{sample}.db.coord-sorted.bam"
    output: "{sample}.db.profile.txt.gz"
    singularity: 'docker://quay.io/arumugamlab/msamtools:MSAM_VERSION_0'
    shell:
        """
        samtools sort -m 20G --threads 4 -n {input} \
            | msamtools filter -b -u -l 80 -p 95 -z 80 --besthit - \
            | msamtools profile --multi=proportional --label={wildcards.sample} -o {output} -
        """
```

The Bioconda container contains the `msamtools` executable, whereas the custom
image also provides the `samtools` executable.

### Container URLs

The exact Bioconda and arumugamlab container image URLs are release-specific
and are listed in the corresponding GitHub release notes. Use those published
URLs rather than guessing the build tag.

## Advanced installation

You can also download the source code and build **msamtools** yourself.

### Required tools

Building from source requires:

1. a C compiler such as `gcc`
2. `make`
3. `gzip`
4. `tar`
5. `bzip2`
6. either `wget` or `curl`

Source builds automatically download and build the required HTSlib release.

### Required libraries

The following development libraries are required:

1. **zlib** development files, e.g. `zlib1g-dev` on Ubuntu
2. **argtable2** development files, e.g. `libargtable2-dev` on Ubuntu
3. pthread support
4. a standard floating-point math library

The latter two are normally provided by the system C environment.

### Building from a release tarball

For normal installation from source, obtain a release tarball and build it
directly:

```bash
wget https://github.com/arumugamlab/msamtools/releases/download/MSAM_VERSION/msamtools-MSAM_VERSION.tar.gz
tar xfz msamtools-MSAM_VERSION.tar.gz
cd msamtools-MSAM_VERSION
./configure
make
make check
```

This builds the `msamtools` executable and runs the deterministic test suite.

To install under the configured installation prefix:

```bash
make install
```

### Building from the Git repository

Developers and contributors can build directly from the Git repository.

#### Getting the source code

Clone the repository:

```bash
git clone https://github.com/arumugamlab/msamtools.git
cd msamtools
```

Repository snapshots can alternatively be downloaded from GitHub, although
cloning with Git is recommended for development work.

#### Bootstrapping the build system

The generated `configure` script is not stored in the Git repository.

Generate the Autotools build files with:

```bash
autoreconf -fi
```

#### Building and testing

Configure, build, and run the test suite with:

```bash
./configure
make
make check
```

The deterministic test suite uses small hand-crafted SAM fixtures and checks
exact behavior of filtering, best-hit and unique-hit selection, profiling,
coverage, alignment summaries, command-line validation, QNAME grouping,
integration between commands, and streaming behavior.

#### Checking a source distribution

Maintainers preparing a source distribution can run:

```bash
make distcheck
```

This creates a source archive, performs an out-of-tree build from that
distribution, and runs the complete test suite against it.

### Building the Docker image

The repository Dockerfile normally builds from a published msamtools release
tarball:

```bash
docker build --build-arg MSAM_VERSION=MSAM_VERSION -t msamtools:MSAM_VERSION .
```

A local release tarball can also be used for testing before publication.
Place `msamtools-<version>.tar.gz` in the Docker build context and run:

```bash
docker build \
    --build-arg MSAM_SOURCE=local \
    --build-arg MSAM_VERSION=<version> \
    -t msamtools:test .
```

The resulting image contains msamtools together with a samtools executable
built to match the HTSlib version used by that msamtools release.

### Validation

Reproducible profile validation using synthetic metagenomic alignments is
described in the [profile validation report](validation/profile_validation.md).

The scripts and input files required for validation are provided under
`validation/`.

The validation results can be reproduced by:

```bash
cd validation
python prepare_genomes.py genomes.tsv
python validate_profiles.py \
    --msamtools ../msamtools \
    --control-community balanced_skew \
    prepared_genomes/validated_genomes.tsv \
    communities/
```

This performs:

1. preparation of the reference genomes;
2. generation of synthetic alignments with known community composition;
3. execution of the msamtools profiling modes;
4. comparison of the resulting profiles with the known truth and generation
   of the validation report and figures.
