# msamtools:  microbiome-related extension to samtools

**msamtools** provides useful functions that are commonly used in microbiome
data analysis, especially when analyzing shotgun metagenomics
or metatranscriptomics data.

## 1. Installation <a name="installation"></a>

### 1.1. System requirements <a name="sys-requirements"></a>

You should be able to use **msamtools** on any flavor of linux and UNIX.
Although I have not tested it myself, it should also work on macOS. A
macOS version is available from bioconda.

**msamtools** has
fixed dependency on **samtools** v1.9, which it will automatically
download and build. samtools has its own requirements even though
I have tuned its configuration to a minimum.

### 1.2. Recommended installation using conda <a name="install-conda"></a>

The easiest way to install **msamtools** and the required dependencies is
via conda from the bioconda channel. Versions 1.0.2 and above are available in bioconda.

If you are already within a conda environment, you can just add it as:
~~~
conda install -c bioconda msamtools
msamtools help
~~~
If you are creating a new environment, then:
~~~
conda create -n msamtools -c conda-forge -c bioconda msamtools
conda activate msamtools
msamtools help
~~~

If you also needed **samtools** for your analysis but it is not available in your path,
you can install them together via conda:
~~~
conda create -n msamtools -c conda-forge -c bioconda msamtools samtools=1.12
conda activate msamtools
msamtools help
samtools help
~~~

If for some reason you cannot install via conda, please check [Advanced installation](#advanced-installation).

### 1.3. Using online docker containers without installing locally<a name="use-docker"></a>

**msamtools** is available as a docker container that can be used e.g. in
snakemake workflows.
There are two possibilities to run **msamtools** using docker containers.

#### 1.3.1. Using a docker container from the BioConda build <a name="use-docker-bioconda"></a>

The first is the **bioconda** docker image corresponding to the bioconda
release. This docker image provides just **msamtools**.
E.g., if you add this line in your snakemake rule
~~~
singularity: 'docker://quay.io/biocontainers/msamtools:1.1.1--h5bf99c6_0'
~~~
you can use this dockerized version of **msamtools** by invoking **snakemake**
as:
~~~
snakemake --use-singularity
~~~

#### 1.3.2. Using a docker container from arumugamlab for msamtools+samtools <a name="use-docker-arumugamlab"></a>

If you need to pipe between **msamtools** and **samtools** (which I do a LOT), then it is
useful to have both **msamtools** and **samtools** in the docker container. Since our conda
release to **bioconda** contains only **msamtools**, we have made a custom container that
contains both **msamtools** and **samtools (v1.9)**. E.g., if you had a bam file that was sorted
by coordinates, which needs to be sorted by name before you can use **msamtools**, you could
have a snakemake rules such as:

~~~
rule profile_sample:
    input: "{sample}.db.coord-sorted.bam"
    output: "{sample}.db.profile.txt.gz"
    singularity: 'docker://quay.io/arumugamlab/msamtools:1.1.1_0'
    shell:
        """
        samtools sort -m 20G --threads 4 -n {input} \\
            | msamtools filter -b -u -l 80 -p 95 -z 80 --besthit - \\
            | msamtools profile --multi=proportional --label={wildcards.sample} -o {output} -
        """
~~~
This will only work with our docker container but not with **bioconda** container.

## 2. Using msamtools <a name="msamtools"></a>

This is the master program that you call with the subprogram options. There
are currently 4 subprograms that you can call as shown below.
~~~

Program: msamtools (Metagenomics-related extension to samtools)
Version: 1.1.1 (using samtools/htslib 1.9)

Usage:   msamtools <command> [options]

Commands:
 -- Filtering
     filter         filter alignments based on alignment statistics

 -- Profiling
     profile        estimate relative abundance profile of reference sequences or genomes in bam file

 -- Coverage
     coverage       estimate per-base or per-sequence read coverage of each reference sequence

 -- Summary
     summary        summarize alignment statistics per read in a table format
~~~

These represent the different analysis of SAM/BAM files you can perform using
**msamtools**.
Section 3 provides example workflows on how to combine them with each other or **samtools**. Sections 4-7 explain the subprograms in detail.

## 3. Example workflows using msamtools<a name="workflows"></a>

Similar to **samtools**, I have designed **msamtools** to work on a stream, avoiding creation of intermediate files.

### 3.1. Alignment and filtering in one step

If your aligner can write to `stdout`, then you can directly pipe the output to **msamtools** and filter on the fly.

#### Task
Align **SAMPLE** (files `SAMPLE.1.fq.gz` and `SAMPLE.2.fq.gz`) to the `bwa-mem2` database in `DB`; retain alignments over `80bp` with `>95%` identity covering `>80%` readlength; and write output to `SAMPLE.DB.filtered.bam`.

#### Command
~~~
bwa-mem2 mem DB SAMPLE.1.fq.gz SAMPLE.2.fq.gz \
   | msamtools filter -S -b -l 80 -p 95 -z 80 > SAMPLE.DB.filtered.bam
~~~

#### Explanation

The command above

* aligns using `bwa-mem2` that generates `SAM` format
* pipes the output to **msamtools**
* asks **msamtools** to
  * read `SAM` format (`-S`)
  * filter alignments that are
    * at least `80bp` long (`-l 80`)
    * at least `95%` identity (`-p 95`)
    * at least `80%` of the read aligned (`-z 80`)
  * write output in `BAM` format (`-b`)

### 3.2. Removing human reads from human metagenomes

Here is an example workflow to filter human reads using `bwa-mem2`.

#### Task
Align **SAMPLE** (fastq files `SAMPLE.1.fq.gz` and `SAMPLE.2.fq.gz`) to the human genome `bwa-mem2` database in `HUMAN_DB`; trust alignments that span at least `30bp`; and write the host-free reads as compressed `fastq` files to `SAMPLE.hostfree.1.fq.gz` and `SAMPLE.hostfree.2.fq.gz`.

>Note: This is our standard workflow for host-sequence removal. Without the extra filtering by `30bp`, we lose significantly more reads that are spuriously flagged as host-derived. Directly going from `bwa-mem2 mem` output to `samtools` is not advisable, as this will remove useful reads from your sample.

#### Command
~~~
bwa-mem2 mem HUMAN_DB SAMPLE.1.fq.gz SAMPLE.2.fq.gz \
  | msamtools filter -S -l 30 --invert --keep_unmapped -bu - \
  | samtools fastq -1 SAMPLE.hostfree.1.fq.gz -2 SAMPLE.hostfree.2.fq.gz -s /dev/null -o /dev/null -c 6 -N -
~~~

#### Explanation

The command above

* aligns using `bwa-mem2` that generates `SAM` format
* pipes the output to **msamtools**
* asks **msamtools** to get reads that are not human by
  * reading `SAM` format (`-S`)
  * filtering alignments that are at least `30bp` long (`-l 30`)
  * negating that and getting alignments that are below `30bp` (`--invert`)
  * while retaining also the unmapped reads (`--keep_unmapped`)
  * writing output in uncompressed `BAM` format (`-bu`)
* then pipes the output to **samtools**
* asks **samtools** to make `fastq` files
  * write compressed forms (`-c 5`)
  * of fastq format (`fastq`)
  * of forward and reverse reads to separate files (`-1` and `-2`)
  * while ignoring unpaired reads (`-s /dev/null -0 /dev/null`)
  * and appending `/1` and `/2` to the reads (`-N`)

### 3.3. Mapping a metagenome sample to a gene database and generating gene profiles

Here is an example workflow one would use after mapping metagenomic reads to IGC.

#### Task

Align **SAMPLE** (fastq files `SAMPLE.1.fq.gz` and `SAMPLE.2.fq.gz`) to the gene catalog `bwa-mem2` database in `GENE_DB`; filter as in Section 3.1 but retain only the highest-scoring hits; and write profile of all genes to `SAMPLE.profile.txt.gz`.

#### Command
~~~
bwa-mem2 mem GENE_DB SAMPLE.1.fq.gz SAMPLE.2.fq.gz \
  | msamtools filter -S -bu -l 80 -p 95 -z 80 --besthit - \
  | msamtools profile --multi=proportional --label=SAMPLE --unit=rel -o SAMPLE.profile.txt.gz -
~~~

#### Explanation

The command above

* aligns using `bwa-mem2` that generates `SAM` format
* then pipes the output to **msamtools filter** * to
  * read `SAM` format (`-S`)
  * filter alignments that are
    * at least `80bp` long (`-l 80`)
    * at least `95%` identity (`-p 95`)
    * at least `80%` of the read aligned (`-z 80`)
  * keep only best-scoring hits per read (`--besthit`)
  * write output in uncompressed `BAM` format (`-bu`)
* then pipes the output to **msamtools profile** * to
  * share multihit reads proportionally among hits (`--multi=proportional`)
  * calculate relative abundance profiles (`--unit=rel`)
  * use sample label `SAMPLE` in the output file (`--label=SAMPLE`)
  * and write compressed output to `SAMPLE.profile.txt.gz` (`-o`)

## 4. msamtools filter <a name="msamtools-filter"></a>

**filter** program provides alignment filtering based on percent identity,
read length, aligned fraction of read length, or combinations thereof.
For example, in mapping metagenomic reads to a database for species-level
annotation, we typically throw out alignments <95% sequence identity.

Here is an example filtering command one would use after mapping metagenomic
reads to the Integrated Gene Catalog (IGC) consisting 9.9 million genes
(Li *et al*, **Nat. biotech** 2014).
~~~
msamtools filter -b -l 80 -p 95 -z 80 --besthit sample1.IGC.bam > sample1.IGC.filtered.bam
~~~
The above command selects all alignments where the read is at least 80bp long,
the percent identity of alignment is >=95%, and at least 80% of the read is
aligned. Additionally, it only keeps the best alignment(s) per read. Multiple
alignments for a read with the same alignment score will be kept.

A full description is given below:
~~~
Usage:
------

msamtools filter [-buhSkv] <bamfile> [--help] [-l <int>] [-p <int>] [--ppt=<int>] [-z <int>] [--rescore] [--besthit] [--uniqhit]

General options:
----------------

These options specify the input/output formats of BAM/SAM files 
(same meaning as in 'samtools view'):
  -b                        output BAM (default: false)
  -u                        uncompressed BAM output (force -b) (default: false)
  -h                        print header for the SAM output (default: false)
  -S                        input is SAM (default: false)
  <bamfile>                 input SAM/BAM file
  --help                    print this help and exit

Specific options:
-----------------

  -l <int>                  min. length of alignment (default: 0)
  -p <int>                  min. sequence identity of alignment, in percentage, integer between 0 and 100; requires NM field to be present (default: 0)
  --ppt=<int>               min/max sequence identity of alignment, in parts per thousand, integer between -1000 and 1000; requires NM field to be present (default: 0)
                            NOTE:
                            -----
                                  When using --ppt, +ve values mean minimum ppt and -ve values mean maximum ppt.
                                  E.g., '--ppt 950' will report alignments with ppt>950,
                                  and '--ppt -950' will report alignments with ppt<=950.
  -z <int>                  min. percent of the query that must be aligned, between 0 and 100 (default: 0)
  -k, --keep_unmapped       report unmapped reads, when filtering using upper-limit thresholds (default: false)
  -v, --invert              invert the effect of the filter (default: false)
                            CAUTION:
                            --------
                                  When using --invert or -v, be precise in what needs to be inverted.
                                  Adding -v gives you the complement of what you get without -v.
                                  Sometimes, this might be counter-intuitive.
                                  E.g., '-l 65 -p 95' will report alignments that are (>65bp AND >95%).
                                        '-l 65 -p 95 -v' will not report (<65bp AND <95%), as one might think.
                                        '-l 65 -p 95 -v' will report NOT(>65bp AND >95%) which is (<65bp OR <95%).
                                        Notice the 'OR' in the final logical operation. This means that
                                        an alignment that fails one condition will still be reported if it
                                        satisfies the other condition.
                                        If you only wanted alignments that are below 95%, then specify '-p 95 -v'
  --rescore                 rescore alignments using MD or NM fields, in that order (default: false)

Special filters:
----------------

The following special filters cannot be combined with -v, but require:
  (1) the alignments to be sorted by name,
  (2) AS field (alignment score) to be present.
You can usually achieve sorting by:
  samtools sort -n -T tmp.sort input.bam  | msamtools -m filter --besthit -
If AS is missing, you can rescore alignments by:
  samtools sort -n -T tmp.sort input.bam | msamtools -m filter --rescore --besthit -

  --besthit                 keep all highest scoring hit(s) per read (default: false)
  --uniqhit                 keep only one highest scoring hit per read, only if it is unique (default: false)
~~~

It is important to note that **filter** command works at the read level.
E.g., using `--besthit` will independently choose the best-scoring alignments
for the forward and reverse reads.

Additionally, when using `--besthit` or `--uniqhit`, the input BAM file
is expected to be sorted by **QNAME**. Each mate should occur as a continuous
group followed by the other mate. Typical output from read-mappers do
produce output that follows this. However, if you had processed the BAM
files yourself and modified things, it could affect output from **filter**.
In such cases, you should sort them again using `samtools sort -n` so they
are sorted by **QNAME** again.

## 5. msamtools profile <a name="msamtools-profile"></a>

**profile** program provides sequence abundance profiling functionality.
By default, relative abundance of each sequence in the BAM file is reported.
However, using the `--genome` option, you can associate database sequences in the
BAM file with larger features, e.g. genomes or MAGs.
Abundance of a sequence/genome is estimated as number of inserts
(mate-pairs/paired-ends) mapped to that sequence/genome, divided by its
length. Reads mapping to multiple sequences/genomes can be shared across
the sequences/genomes in three different ways (please see below).
Finally, abundance is estimated in one of four units:
abundance (ab), relative abundance (rel),
fragments per kilobase of sequence per million reads (fpkm),
or transcripts per million (tpm). As you probably understand, *tpm* and
*fpkm* are probably not suitable for profiling genomes, but do not let me
stop you!

>**WARNING: The profiler expects that BAM files are sorted by name so that
it can keep track of reads that map to multiple locations. Please ensure
that your BAM files are sorted that way. Profiler does not check this, so
can give you erroneous results when you pass coordinate-sorted BAM files.**


>**NOTE:** From **v1.0.0**, the default output is a gzipped text file. Therefore,
argument `--gzip` or `-z` will throw an error.

We highly recommend that you filter the alignments before sending to the
**profile** program, as it considers each alignment to be important (it
does not look at alignment quality, for example).

Here is an example profiling command one would use after mapping metagenomic
reads to IGC.
~~~
msamtools profile --multi=proportional --label=sample1 --unit=rel -o sample1.IGC.profile.txt.gz sample1.IGC.filtered.bam
~~~
The above command estimates relative abundance of IGC genes after sharing
multi-mapper reads proportionately between the genes (see below).

In the spirit of **samtools** programs, **msamtools** programs can also
stream between each other. Therefore, a single command to **filter** and **profile**
would look like:
~~~
msamtools filter -b -u -l 80 -p 95 -z 80 --besthit sample1.IGC.bam \
  | msamtools profile --multi=proportional --label=sample1 --unit=rel -o sample1.IGC.profile.txt.gz -
~~~

or for mapping to scaffolds that are grouped in metagenome-assembled MAGs using:
~~~
msamtools filter -b -u -l 80 -p 95 -z 80 --besthit sample1.myMAGs.bam \
  | msamtools profile --multi=proportional --label=sample1 --unit=rel --genome myMAGs.genome.def -o sample1.myMAGs.profile.txt.gz -
~~~

### 5.1. Profiling genomes or MAGs <a name="profiling-genomes"></a>

Starting from **v1.0.0**, **profile** program supports profiling of genomes defined by a set of
sequences. This requires a tab-delimited definition file of the following format:
~~~
MAG_1	Contig_1
MAG_1	Contig_2
MAG_1	Contig_3
MAG_1	Contig_4
MAG_2	Contig_18
MAG_2	Contig_27
MAG_2	Contig_32
MAG_3	Contig_23
MAG_3	Contig_24
MAG_3	Contig_35
MAG_3	Contig_48
~~~

If this information is stored in a file called `myMAGs.genome.def`, then you can
run the profiler as follows to get profiles at the genome level.
~~~
msamtools filter -b -u -l 80 -p 95 -z 80 --besthit sample1.myMAGs.bam \
  | msamtools profile --multi=proportional --label=sample1 --unit=rel --genome myMAGs.genome.def -o sample1.myMAGs.profile.txt.gz -
~~~

### 5.2. Units of abundance <a name="abundance-units"></a>

By default, **profile** command will generate relative abundances that sum to `1`
across the sequences in the BAM file. Four options to measure
the abundance are available:
*  **ab** - number of inserts mapped to the sequence, normalized by sequence length
*  **rel** - relative abundance, which is **ab** normalized by sum across all sequences
* **fpkm** - fragments per kilobase of sequence per million reads
* **tpm** - transcripts per million

An optional `--nolen` flag turns off sequence length normalization for **ab** and **rel**.
When combining `--unit=ab` and ` --nolen`, you get the raw number of inserts mapped
to each sequence, and summing them up will match the total number of inserts in
the BAM file (or what was passed via `--total`).

### 5.3. Keeping track of unmapped reads <a name="track-unmapped"></a>

By default, **profile** command will generate relative abundances that sum to `1`
across the sequences in the BAM file. In metagenomic data, sometimes we need
to identify the fraction of the reads that were not mapped to our database
and only assign the remaining fraction to the sequences in the BAM file.
Consider the following examples. Please note that in the examples below we do
not address issues inherent to relative abundance, which is beyond the scope of
a profiler, but merely provide different ways to estimate it.

**Sample 1** has 500 cells, with abundance given in column 2.
If we sequenced 1000 reads from this samples, we expect to see
the number of reads in column 3 (we ignore the
difference in genome size for this simple example).
The relative abundance of the taxa estimated when **excluding**
the unmapped reads is given in column 4.
The relative abundance of the taxa estimated when **including**
the unmapped reads is given in column 5.

**Sample 1:**

| Taxon            | Abundance | Reads | -Unknown | +Unknown |
|------------------|----------:|------:|---------:|---------:|
| Bacteroides      |       100 |  200  |    25%   |   20%    |
| Prevotella       |       200 |  400  |    50%   |   40%    |
| Faecalibacterium |       100 |  200  |    25%   |   20%    |
| Unknown          |       100 |  200  |    ---   |   20%    |

Here's another sample, **Sample 2**, with different microbial composition.

**Sample 2:**

| Taxon            | Abundance | Reads | -Unknown | +Unknown |
|------------------|----------:|------:|---------:|---------:|
| Bacteroides      |       100 |  200  |    50%   |   20%    |
| Faecalibacterium |       100 |  200  |    50%   |   20%    |
| Unknown          |       300 |  600  |    ---   |   60%    |

Even though *Bacteroides* and *Faecalibacterium* had the same abundance
in the two samples, ignoring the **Unknown** fraction leads to a different
estimation of their relative abundances. In these situations, it is useful
to keep track of the **unmapped** reads from the metagenome.

If you would like to include the **Unknown** fraction in the profile, you should
tell the profiler how many reads were originally sequenced. Then the profiler will estimate
the unmapped reads based on how many reads are in the bam/sam file, and then use
it in the profiling stage.

~~~
# Get number of entries in the fwd fastq file = number of inserts
lines=$(zcat sample1.1.fq.gz | wc -l)
entries=$(expr \$lines / 4)   # There are 4 lines per fastq entry

# Use total reads in profiler
msamtools filter -b -u -l 80 -p 95 -z 80 --besthit sample1.IGC.bam \
  | msamtools profile --multi=proportional --label=sample1 -o sample1.IGC.profile.txt.gz --total=$entries -
~~~

### 5.4. Avoiding extremely-low-abundant features <a name="profile-mincount"></a>

When just a handful of reads map to a feature (genome or contig or gene) in the database,
it is not immediately clear if it is just really low abundance or if it was a spurious
mapping. While it might be real rare features, it is sometimes preferable to only
take a feature forward to downstream analysis when a reasonable number of reads map to it.
From **v1.1.0**, you can use `--mincount` to specify the minimum number of reads that a
feature should attract
for it to be considered **detected** - meaning **expressed** in metatranscriptomic data or
**present** in metagenomic data. The specific threshold should be based on the sequencing
depth of the sample. While the default behavior is to not apply this filter, we
recommend to use `10` for metagenomes or metatranscriptomes with `>10M` paired-end
reads.

### 5.5. Useful information in the output file <a name="profile-output"></a>

The header section of the output file includes a few lines of comment that
are hopefully useful. Here is an example:

~~~
# msamtools version 1.1.0
# Command: msamtools profile --label test --unit rel --multi prop --total 3519692 -o test.profile.txt test.bam
#   Total inserts: 3519692
#  Mapped inserts: 334063
# Mapped fraction: 0.0949
# Estimated seq. length for 'Unknown': 6234bp
001-02
Unknown 0.809179
001-02_NODE_90_length_5676_cov_189.052854   0.0109423
001-04_NODE_36_length_4618_cov_68.000866    0
...
~~~

The commented lines are self-explanatory, and could be useful in getting
quick summary of the profiling process. Since length-normalization is not
turned off, the average sequence length for the entire database is used
as a proxy sequence length for the **Unknown** fraction.

The first line includes the name of the sample provided via `--label`.
This is for conveniently combining output from multiple files. A script that
combines output does not need external information to create a table with the
right sample name in a row/column.

A full description is given below:
~~~
Usage:
------

msamtools profile [-S] <bamfile> [--help] -o <file> --label=<string> [--genome=<string>] [--total=<int>] [--mincount=<int>] [--unit=<string>] [--nolen] [--multi=<string>]

General options:
----------------

These options specify the input/output formats of BAM/SAM files 
(same meaning as in 'samtools view'):
  -S                        input is SAM (default: false)
  <bamfile>                 input SAM/BAM file
  --help                    print this help and exit

Specific options:
-----------------

  -o <file>                 name of output file (required)
  --label=<string>          label to use for the profile; typically the sample id (required)
  --genome=<string>         tab-delimited genome definition file - 'genome-id<tab>seq-id' (default: none)
  --total=<int>             number of high-quality inserts (mate-pairs/paired-ends) that were input to the aligner (default: 0)
  --mincount=<int>          minimum number of inserts mapped to a feature, below which the feature is counted as absent (default: 0)
  --unit=<string>           unit of abundance to report {ab | rel | fpkm | tpm} (default: rel)
  --nolen                   do not normalize the abundance (only relevant for ab or rel) for sequence length (default: normalize)
  --multi=<string>          how to deal with multi-mappers {all | equal | proportional} (default: proportional)

Description
-----------

Produces an abundance profile of all reference sequences in a BAM file
based on the number of read-pairs (inserts) mapping to each reference sequence.
It can work with genome-scale reference sequences while mapping to a database 
of sequenced genomes, but can also work with gene-scale sequences such as in the
Integrated Gene Catalog from human gut microbiome (Li et al, Nat biotech 2014).

In the output file, each sequence in the BAM file gets a line with its abundance.
They are presented in the order in which they appear in the BAM header. <label>
is used as the first line, so that reading or 'joining' these files is easier.

--total option:      In metagenomics, an unmapped read could still be a valid
                     sequence, just missing in the database being mapped against.
                     This is the purpose of the '--total' option to track the
                     fraction of 'unknown' entities in the sample. If --total
                     is ignored or specified as --total=0, then tracking the 
                     'unknown' fraction is disabled. However, if the total 
                     sequenced inserts were given, then there will be a new
                     feature added to denote the 'unknown' fraction.
Units of abundance:  Currently four different units are available.
                          'ab': raw insert-count abundance
                         'rel': relative abundance (default)
                        'fpkm': fragments per kilobase of sequence per million reads
                         'tpm': transcripts per million
                     If number of reads input to the aligner is given via --total,
                     fpkm and tpm will behave differently than in RNAseq data,
                     as there is now a new entity called 'unknown'.
Alignment filtering: 'profile' expects that every alignment listed is considered 
                     valid. For example, if one needs to filter alignments 
                     based on alignment length, read length, alignment percent
                     identity, etc, this should have been done prior to 
                     'profile'. Please see 'filter' for such filtering.
Multi-mapper reads:  Reads mapping to multiple references need to be considered
                     carefully, as spurious mappings of promiscuous regions or
                     short homology could lead to incorrect abundances of 
                     sequences. 'profile' offers three options for this purpose.
                     If an insert maps to N references at the same time:
                   'all': each reference gets 1 insert added.
                 'equal': each reference gets 1/N insert added.
          'proportional': each reference gets a fraction proportional to its 
                          reference-sequence-length-normalized relative 
                          abundance estimated only based on uniquely
                          mapped reads.
~~~

## 6. msamtools coverage <a name="msamtools-coverage"></a>

**coverage** program estimates per-position or fractional coverage of each sequence in
the BAM file.

### 6.1. Per-position coverage of all sequences in BAM file <a name="pos-coverage"></a>

The per-position coverage output file is in the format of old Sanger quality files
with fasta headers and space-delimited numbers. As this can build up into
quite a large file, using the `-x` option will not print coverage for
sequences that did not have a single read mapped to them. Since their coverage
is essentially zero in each position, printing their coverage is just a
waste of space.

Here is an example per-position coverage command.
~~~
msamtools coverage -x -z -o sample1.coverage.txt.gz sample1.IGC.filtered.bam
~~~

### 6.2. Fractional coverage of each sequence in BAM file <a name="frac-coverage"></a>

Sometimes it is useful to see which sequence from the BAM file has been observed in
the sample. And if yes, it is nice to know what fraction of the sequence has been
covered with alignments in the BAM file. For this one can use the `--summary`
option, which outputs fractional coverage and sequencing-coverage of each sequence.

Here is an example fractional coverage command.
~~~
msamtools coverage -z --summary -o sample1.coverage.summary.txt.gz sample1.IGC.filtered.bam
~~~

And here is an example output:
~~~
cluster_001_consensus_length_3171293	0.05464883	0.25
cluster_002_consensus_length_2788722	0.99955930	10.79
cluster_003_consensus_length_6395848	0.99998921	38.10
cluster_004_consensus_length_2025181	0.99947906	31.14
cluster_005_consensus_length_3532514	0.99987346	70.04
~~~
First column names the sequence, 2nd column reports the fraction of that sequence that is covered
and the 3rd column gives sequencing-coverage. Apparently, the 5th genome has 70X coverage in
that sample!

A full description is given below:
~~~
Usage:
------

msamtools coverage [-Sxz] <bamfile> [--help] -o <file> [--summary] [-w <int>]

General options:
----------------

These options specify the input/output formats of BAM/SAM files 
(same meaning as in 'samtools view'):
  -S                        input is SAM (default: false)
  <bamfile>                 input SAM/BAM file
  --help                    print this help and exit

Specific options:
-----------------

  -o <file>                 name of output file (required)
  --summary                 do not report per-position coverage but report fraction of sequence covered (default: false)
  -x, --skipuncovered       do not report coverage for sequences without aligned reads (default: false)
  -w, --wordsize=<int>      number of words (coverage values) per line (default: 17)
  -z, --gzip                compress output file using gzip (default: true)

Description:
------------

Produces per-position sequence coverage information for all reference sequences
in the BAM file. Output is similar to old-style quality files from the Sanger 
sequencing era, with a fasta-style header followed by lines of space-delimited 
numbers.

For large datasets, option '-x' comes in handy when only a small fraction of 
reference sequences are covered.

If using '-z', output file does NOT automatically get '.gz' extension. This is 
up to the user to specify the correct full output file name.
~~~

## 7. msamtools summary <a name="msamtools-summary"></a>

**summary** program summarizes alignments given in the BAM file. It can
also provide distributions of certain features across alignments.

Here is an example summary command that reports all alignments:
```console
mani@host:~/src/git/msamtools> ./msamtools summary input.bam | head
ERR688505.1_FCD1R4CACXX:7:1101:1787:2090#CAGCGGCG       93      MH0153_GL0123525        93      93      100.0
ERR688505.1_FCD1R4CACXX:7:1101:1787:2090#CAGCGGCG       73      MH0153_GL0123525        73      73      100.0
ERR688505.2_FCD1R4CACXX:7:1101:1892:2123#CAGCGGCG       93      MH0455_GL0055430        93      93      100.0
ERR688505.2_FCD1R4CACXX:7:1101:1892:2123#CAGCGGCG       89      MH0060_GL0037700        89      60      67.4
ERR688505.3_FCD1R4CACXX:7:1101:1752:2179#CAGCGGCN       92      O2.UC22-2_GL0054026     92      33      35.9
ERR688505.3_FCD1R4CACXX:7:1101:1752:2179#CAGCGGCN       92      O2.UC22-2_GL0054026     92      92      100.0
ERR688505.4_FCD1R4CACXX:7:1101:1788:2199#CAGCGGCG       93      MH0204_GL0114410        93      91      97.8
ERR688505.4_FCD1R4CACXX:7:1101:1788:2199#CAGCGGCG       62      MH0204_GL0114410        62      62      100.0
ERR688505.5_FCD1R4CACXX:7:1101:1765:2211#CAGCGGCN       91      MH0188_GL0130879        91      91      100.0
ERR688505.5_FCD1R4CACXX:7:1101:1765:2211#CAGCGGCN       78      MH0188_GL0130879        78      78      100.0
```

Here is another example summary command that reports the distribution
of unmapped bases:
```console
mani@host:~/src/git/msamtools> ./msamtools summary --stats=unmapped input.bam
0       50033
1       24425
2       12495
3       6548
4       3426
5       1913
6       1241
7       828
8       609
9       499
10      481
11      432
12      405
13      313
14      327
15      339
16      337
17      318
18      301
19      304
20      285
21      256
22      261
23      279
24      257
25      272
26      269
27      259
28      265
29      243
30      265
31      254
32      241
33      261
34      248
35      220
36      226
37      229
38      258
39      241
40      223
41      239
42      236
43      239
44      236
45      257
46      226
47      229
48      246
49      250
50      205
51      238
52      234
53      200
54      218
55      226
56      227
57      208
58      212
59      195
60      166
61      158
62      167
63      174
64      12
65      11
66      9
67      4
68      6
69      2
70      5
71      3
72      2
73      5
74      1
```
It shows that 50033 reads had 0 unmapped bases (meaning full mapping), 24425 reads had 1 unmapped base, etc.

A full description is given below:
~~~
Usage:
------

msamtools summary [-Sc] <bamfile> [--help] [-e <num>] [--stats=<string>]

General options:
----------------

These options specify the input/output formats of BAM/SAM files 
(same meaning as in 'samtools view'):
  -S                        input is SAM (default: false)
  <bamfile>                 input SAM/BAM file
  --help                    print this help and exit

Specific options:
-----------------

  -e, --edge=<num>          ignore alignment if reads map to <num> bases at the edge of target sequence (default: 0)
  -c, --count               count number of unique inserts in BAM file (default: false)
  --stats=<string>          {mapped|unmapped|edit|score} only report readcount distribution for specified stats, not read-level stats (default: none)

Description
-----------
Prints summary of alignments in the given BAM/SAM file. By default, it prints
a summary line per alignment entry in the file. The summary is a tab-delimited
line with the following fields:
	qname,aligned_qlen,target_name,glocal_align_len,matches,percent_identity
glocal_align_len includes the unaligned qlen mimicing a global alignment 
in the query and local alignment in target, thus glocal.

With --stats option, summary is consolidated as distribution of read counts
for a given measure. 
   --stats=mapped   - distribution for number of mapped query bases
   --stats=unmapped - distribution for number of unmapped query bases
   --stats=edit     - distribution for edit distances
   --stats=score    - distribution for score=match-edit
~~~

## 8. Advanced installation <a name="advanced-installation"></a>

You can also download the source code and build it yourself.

### 8.1. Required tools <a name="required-tools"></a>

While building **msamtools**, you will need some standard tools that are
most likely installed in your system by default. I will still list them here
anyway to be sure:

 1. gcc
 2. gzip
 3. tar
 4. wget

If any of these is missing in your system, or cannot be found in your
application path, please fix that first.

### 8.2. Required libraries <a name="required-libraries"></a>

The following libraries are required to build **msamtools** from source:

 1. **zlib** development version (e.g., zlib1g-dev in ubuntu)
 2. **argtable2** development version (e.g., libargtable2-dev in ubuntu)

Please make sure that these are installed in your system before trying to
build.

### 8.3. For normal users <a name="normal-users"></a>

If you are a normal user, then the easiest way is to obtain the package file
and build the program right away. The following commands were written when
version 1.1.0 was the latest, so please update the version number in the
commands below.

**Note:** Newer C compilers from gcc use `-std=gnu99` by default, which I had
not tested on version 0.9 as my gcc version is quite outdated with `-std=gnu89` as default.
This leads to version 0.9 not compiling when running `make` with new compilers. The
current fix for using the release tarball for version 0.9 is to tell the compiler which
standard to use, using `CFLAGS="-std=gnu89"`. This extra option was a
temporary fix only, and is not needed from version 0.9.6 as the code has been
upgraded to be compatible with `-std=gnu99`.
(Thanks [Russel88](https://github.com/Russel88) for reporting this).

```console
wget https://github.com/arumugamlab/msamtools/releases/download/1.1.0/msamtools-1.1.0.tar.gz
tar xfz msamtools-1.1.0.tar.gz
cd msamtools-1.1.0
./configure
make
```

This should create `msamtools` executable.

### 8.4. For advanced users <a name="advanced-users"></a>

If you are an advanced user who would like to contribute to the code base
or if you just like to do things the hard way, you can check out the source
code and build the program in a series of steps involving `autoconf` and
`automake`. If these names confuse you or scare you, then please follow the
instructions for [normal users](#normal-users).

#### 8.4.1. Getting the source code <a name="source-code"></a>

You can get **msamtools** code from github at
<https://github.com/arumugamlab/msamtools>.
You can either `git clone` it or download the ZIP file and extract the
package.

##### 8.4.1.1. Cloning the git repository <a name="git-clone"></a>

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

##### 8.4.1.2. Downloading the ZIP file from github <a name="git-zip"></a>

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

#### 8.4.2. Running autoconf and automake <a name="automake"></a>

You can check the contents of the repository in the package directory.
```console
$ ls
configure.ac  make_readme.sh  mMatrix.c        msam_helper.c   tests
deps          mBamVector.c    mMatrix.h        msam_profile.c  versions.txt
Dockerfile    mBamVector.h    msam_coverage.c  msam_summary.c  zoeTools.c
LICENSE       mCommon.c       msam_filter.c    msamtools.c     zoeTools.h
Makefile.am   mCommon.h       msam.h           README.md
```

You will note that the `configure` script does not exist in the package.
This is because you need to generate the `configure` script using
`aclocal`, `autoconf` and `automake`.
```
aclocal
autoconf
mkdir build-aux
automake --add-missing
```

#### 8.4.3. Building the program <a name="build"></a>

You can then build msamtools as follows:
```console
./configure
make
```

