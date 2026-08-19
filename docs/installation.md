# Installation

## System requirements <a name="sys-requirements"></a>

You should be able to use **msamtools** on any flavor of linux and UNIX.
Although I have not tested it myself, it should also work on macOS. A
macOS version is available from bioconda.

**msamtools** has
fixed dependency on **samtools** v1.9, which it will automatically
download and build. samtools has its own requirements even though
I have tuned its configuration to a minimum.

## Recommended installation using conda <a name="install-conda"></a>

The easiest way to install **msamtools** and the required dependencies is
via conda from the bioconda channel. Versions 1.0.2 and above are available in bioconda.

If you are already within a conda environment, you can just add it as:
```bash
conda install -c bioconda msamtools
msamtools help
```
If you are creating a new environment, then:
```bash
conda create -n msamtools -c conda-forge -c bioconda msamtools
conda activate msamtools
msamtools help
```

If you also needed **samtools** for your analysis but it is not available in your path,
you can install them together via conda:
```bash
conda create -n msamtools -c conda-forge -c bioconda msamtools samtools=1.12
conda activate msamtools
msamtools help
samtools help
```

If for some reason you cannot install via conda, please check [Advanced installation](#advanced-installation).

## Using online docker containers without installing locally<a name="use-docker"></a>

**msamtools** is available as a docker container that can be used e.g. in
`Snakemake` workflows.
There are two possibilities to run **msamtools** using docker containers.

### Using a docker container from the BioConda build <a name="use-docker-bioconda"></a>

The first is the **bioconda** docker image corresponding to the bioconda
release. This docker image provides just **msamtools**.
E.g., if you add this line in your snakemake rule
```snakemake
singularity: 'docker://quay.io/biocontainers/msamtools:1.1.3--h577a1d6_1'
```
you can use this dockerized version of **msamtools** by invoking **snakemake**
as:
```bash
snakemake --use-singularity
```

### Using a docker container from arumugamlab for msamtools+samtools <a name="use-docker-arumugamlab"></a>

If you need to pipe between **msamtools** and **samtools** (which I do a LOT), then it is
useful to have both **msamtools** and **samtools** in the docker container. Since our conda
release to **bioconda** contains only **msamtools**, we have made a custom container that
contains both **msamtools** and **samtools (v1.9)**. E.g., if you had a bam file that was sorted
by coordinates, which needs to be sorted by name before you can use **msamtools**, you could
have a snakemake rules such as:

```snakemake
rule profile_sample:
    input: "{sample}.db.coord-sorted.bam"
    output: "{sample}.db.profile.txt.gz"
    singularity: 'docker://quay.io/arumugamlab/msamtools:1.1.3_0'
    shell:
        """
        samtools sort -m 20G --threads 4 -n {input} \\
            | msamtools filter -b -u -l 80 -p 95 -z 80 --besthit - \\
            | msamtools profile --multi=proportional --label={wildcards.sample} -o {output} -
        """
```
This will only work with our docker container but not with **bioconda** container.

# Advanced installation

You can also download the source code and build it yourself.

## Required tools <a name="required-tools"></a>

While building **msamtools**, you will need some standard tools that are
most likely installed in your system by default. I will still list them here
anyway to be sure:

 1. gcc
 2. gzip
 3. tar
 4. wget
 5. bzip2 (for installing HTSlib)

If any of these is missing in your system, or cannot be found in your
application path, please fix that first.

## Required libraries <a name="required-libraries"></a>

The following libraries are required to build **msamtools** from source:

 1. **zlib** development version (e.g., zlib1g-dev in ubuntu)
 2. **argtable2** development version (e.g., libargtable2-dev in ubuntu)

Please make sure that these are installed in your system before trying to
build.

## For normal users <a name="normal-users"></a>

If you are a normal user, then the easiest way is to obtain a release
tarball and build the program directly from it.

**Note:** Newer C compilers from gcc use `-std=gnu99` by default, which I had
not tested on version 0.9 as my gcc version is quite outdated with `-std=gnu89` as default.
This leads to version 0.9 not compiling when running `make` with new compilers. The
current fix for using the release tarball for version 0.9 is to tell the compiler which
standard to use, using `CFLAGS="-std=gnu89"`. This extra option was a
temporary fix only, and is not needed from version 0.9.6 as the code has been
upgraded to be compatible with `-std=gnu99`.
(Thanks [Russel88](https://github.com/Russel88) for reporting this).

```bash
wget https://github.com/arumugamlab/msamtools/releases/download/MSAM_VERSION/msamtools-MSAM_VERSION.tar.gz
tar xfz msamtools-MSAM_VERSION.tar.gz
cd msamtools-MSAM_VERSION
./configure
make
make check
```

This builds the `msamtools` executable and runs the fast deterministic test
suite. To install the program under the configured installation prefix, run:
```bash
make install
```

## For advanced users <a name="advanced-users"></a>

If you are an advanced user who would like to contribute to the code base,
you can check out the source code, bootstrap the build system using GNU
Autotools, and build and test the program. If these tools are unfamiliar,
please follow the instructions for [normal users](#normal-users).

### Getting the source code <a name="source-code"></a>

You can get **msamtools** code from github at
<https://github.com/arumugamlab/msamtools>.
You can either `git clone` it or download the ZIP file and extract the
package.

#### Cloning the git repository <a name="git-clone"></a>

You can get a clone of the repository if you wish to keep it up-to-date
independent of our releases.

```console
$ git clone https://github.com/arumugamlab/msamtools.git
Cloning into 'msamtools'...
remote: Enumerating objects: 285, done.
remote: Counting objects: 100% (285/285), done.
remote: Compressing objects: 100% (181/181), done.
remote: Total 285 (delta 167), reused 215 (delta 101), pack-reused 0
Receiving objects: 100% (285/285), 130.93 KiB | 0 bytes/s, done.
Resolving deltas: 100% (167/167), done.
$ cd msamtools
```

You can check the contents of the repository in *msamtools* directory.

#### Downloading the ZIP file from github <a name="git-zip"></a>

You can download the repository snapshot as on the day of download by:
```console
$ wget https://github.com/arumugamlab/msamtools/archive/master.zip
--2021-11-17 12:24:24--  https://github.com/arumugamlab/msamtools/archive/master.zip
Resolving github.com (github.com)... 140.82.121.4
Connecting to github.com (github.com)|140.82.121.4|:443... connected.
HTTP request sent, awaiting response... 302 Found
Location: https://codeload.github.com/arumugamlab/msamtools/zip/master [following]
--2021-11-17 12:24:25--  https://codeload.github.com/arumugamlab/msamtools/zip/master
Resolving codeload.github.com (codeload.github.com)... 140.82.121.10
Connecting to codeload.github.com (codeload.github.com)|140.82.121.10|:443... connected.
HTTP request sent, awaiting response... 200 OK
Length: unspecified [application/zip]
Saving to: ‘master.zip’

     0K .......... .......... .......... .......... .......... 1.68M
    50K .......... ........                                    67.2M=0.03s

2021-11-17 12:24:25 (2.28 MB/s) - ‘master.zip’ saved [70091]

$ unzip master.zip
$ cd msamtools-master
```

### Bootstrapping the build system <a name="automake"></a>

You can check the contents of the repository in the package directory.
```console
$ ls
configure.ac  make_readme.sh  mMatrix.c        msam_helper.c   tests
deps          mBamVector.c    mMatrix.h        msam_profile.c  versions.txt
Dockerfile    mBamVector.h    msam_coverage.c  msam_summary.c  zoeTools.c
LICENSE       mCommon.c       msam_filter.c    msamtools.c     zoeTools.h
Makefile.am   mCommon.h       msam.h           README.md
```

You will note that the `configure` script does not exist in the repository.
Generate the Autotools build files with:
```bash
autoreconf -fi
```

### Building the program <a name="build"></a>

You can then build and test msamtools as follows:
```bash
./configure
make
make check
```

The fast deterministic test suite uses small hand-crafted SAM fixtures and
checks exact behavior of filtering, best-hit and unique-hit selection,
profiling, coverage, alignment summaries, command-line validation, QNAME
grouping, integration between commands, and streaming behavior.

### Checking a source distribution <a name="distcheck"></a>

Maintainers preparing a source distribution can run:
```bash
make distcheck
```

This creates a source archive, performs an out-of-tree build from that
distribution, and runs the complete test suite against it.

### Validation <a name="validation"></a>

Reproducible profile validation using synthetic metagenomic alignments is
described in the [profile validation report](validation/profile_validation.md).

The scripts and input files needed for these steps are provided under
`validation/`. See the profile validation report for the
parameters and validation results.

The validation results can be reproduced by:

```bash
cd validation
python prepare_genomes.py genomes.tsv
python validation/validate_profiles.py --msamtools ../msamtools --control-community balanced_skew prepared_genomes/validated_genomes.tsv communities/
```

This does the following:
1. preparing the reference genomes;
2. generating synthetic alignments with known community composition;
3. running the msamtools profiling modes;
4. comparing the resulting profiles with the known truth and generating the
   validation report and figures.
