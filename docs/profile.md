# msamtools profile

The `profile` command provides sequence abundance profiling functionality.
By default, relative abundance of each sequence in the BAM file is reported.
However, using the `--genome` option, you can associate database sequences in
the BAM file with larger features, e.g. genomes or MAGs.

Abundance of a sequence/genome is estimated from the number of inserts
(mate-pairs/paired-ends) mapped to that sequence/genome, with sequence-length
normalization applied by default. Inserts mapping to multiple sequences/genomes
can be handled in four different ways (see below).

Finally, abundance is estimated in one of four units:
abundance (`ab`), relative abundance (`rel`),
fragments per kilobase of sequence per million reads (`fpkm`),
or transcripts per million (`tpm`).

Here is an example profiling command one would use after mapping metagenomic
reads to IGC.

```bash
msamtools profile --multi=proportional --label=sample1 --unit=rel -o sample1.IGC.profile.txt.gz sample1.IGC.filtered.bam
```

The command estimates relative abundance of IGC genes after sharing
multi-mapper inserts proportionally between the genes (see below).

See [Example workflows](examples.md) for examples of streaming `filter` and
`profile` together.

## Important usage notes

- `profile` requires contiguous QNAME groups to identify inserts and
  multi-mapping correctly. See
  [Input requirements](input-requirements.md#qname-grouping).

- From v1.0.0, `profile` output is always gzip-compressed. The `--gzip` and
  `-z` arguments are therefore not accepted.

- From v1.2.0, pandas-style output is the default. Use `--no-pandas` to obtain
  the previous default output format.

- We highly recommend filtering alignments using `msamtools filter` before sending them to `profile`,
  because each input alignment is treated as valid and `profile` does not, for
  example, apply alignment-quality filtering itself.

## Profiling genomes or MAGs

Starting from **v1.0.0**, the `profile` command supports profiling of genomes
defined by sets of sequences. This requires a tab-delimited definition file of
the following format:

```text
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
```

If this information is stored in a file called `myMAGs.genome.def`, then you
can run the profiler as follows to obtain profiles at the genome level:

```bash
msamtools filter -b -u -l 80 -p 95 -z 80 --besthit sample1.myMAGs.bam \
  | msamtools profile --multi=proportional --label=sample1 --unit=rel --genome myMAGs.genome.def -o sample1.myMAGs.profile.txt.gz -
```

## Units of abundance

By default, the `profile` command generates relative abundances that sum to
`1` across the reported features.

Four abundance units are available:

* **ab** - insert-count abundance, normalized by sequence length by default
* **rel** - relative abundance, obtained by normalizing **ab** across all features
* **fpkm** - fragments per kilobase of sequence per million reads
* **tpm** - transcripts per million

The optional `--nolen` flag turns off sequence-length normalization for
**ab** and **rel**.

When combining `--unit=ab` and `--nolen`, the reported feature values are
insert-equivalents assigned to each feature before sequence-length
normalization. Their sum is reported in the profile header as **Effective
insert-equivalents**.

How this quantity relates to the number of physical mapped inserts depends on
the selected multi-mapper handling mode:

* With `--multi=ignore`, only uniquely mapped inserts contribute to the profile.
* With `--multi=equal`, each successfully assigned multi-mapping insert
  contributes fractions that sum to one insert-equivalent across its matched
  features.
* With `--multi=proportional`, each successfully assigned multi-mapping insert
  likewise contributes fractions that sum to one insert-equivalent, but the
  fractions are assigned according to the iteratively estimated feature
  abundances.
* With `--multi=all`, each matched feature receives one full
  insert-equivalent. A single physical insert can therefore contribute more
  than one insert-equivalent, and the summed profile abundance can exceed both
  the number of mapped inserts and the total number of input inserts.

If `--mincount` removes low-abundance features, their assigned contributions
are removed from the retained profile and counted as **Purged
insert-equivalents**. The remaining feature values still sum exactly to the
reported **Effective insert-equivalents**.

> **NOTE:** FPKM and TPM are conventionally used for gene- or transcript-level
> abundance. For genome/MAG profiling, consider whether `ab` or `rel` is more
> appropriate for the intended downstream analysis.

## Keeping track of unmapped reads

By default, the `profile` command generates relative abundances that sum to
`1` across the features represented in the BAM file. In metagenomic data,
however, it can be useful to retain information about the fraction of the
sequenced inserts that did not map to the reference database.

Consider the following examples. Please note that these examples do not
address issues inherent to relative abundance, which are beyond the scope of
the profiler, but merely illustrate different ways to estimate it.

**Sample 1** has 500 cells, with abundance given in column 2.
If we sequenced 1000 reads from this sample, we expect to see the number of
reads in column 3 (ignoring differences in genome size for this simple
example). The relative abundance estimated when **excluding** the unmapped
reads is given in column 4, while the estimate when **including** them is
given in column 5.

**Sample 1:**

| Taxon            | Abundance | Reads | -Unknown | +Unknown |
|------------------|----------:|------:|---------:|---------:|
| Bacteroides      |       100 |   200 |      25% |      20% |
| Prevotella       |       200 |   400 |      50% |      40% |
| Faecalibacterium |       100 |   200 |      25% |      20% |
| Unknown          |       100 |   200 |      --- |      20% |

Here is another sample, **Sample 2**, with a different microbial composition.

**Sample 2:**

| Taxon            | Abundance | Reads | -Unknown | +Unknown |
|------------------|----------:|------:|---------:|---------:|
| Bacteroides      |       100 |   200 |      50% |      20% |
| Faecalibacterium |       100 |   200 |      50% |      20% |
| Unknown          |       300 |   600 |      --- |      60% |

Even though *Bacteroides* and *Faecalibacterium* had the same abundance in
the two samples, ignoring the **Unknown** fraction leads to different relative
abundance estimates. In these situations, it is useful to keep track of the
unmapped inserts from the metagenome.

To include the **Unknown** fraction in the profile, provide the total number of
high-quality inserts that were input to the aligner using `--total`.
When supplied, `--total` must be a positive integer.

For example:

```bash
# Get number of entries in the forward FASTQ file = number of inserts
lines=$(zcat sample1.1.fq.gz | wc -l)
entries=$(expr $lines / 4)   # There are 4 lines per FASTQ entry

# Use total inserts in profiler
msamtools filter -b -u -l 80 -p 95 -z 80 --besthit sample1.IGC.bam \
  | msamtools profile --multi=proportional --label=sample1 -o sample1.IGC.profile.txt.gz --total=$entries -
```

If `--total` is omitted, tracking of the **Unknown** fraction is disabled.

## Avoiding extremely-low-abundant features

When only a small number of inserts are assigned to a feature (genome,
contig, or gene), it may not be clear whether the feature is genuinely present
at very low abundance or whether the assignments result from spurious mapping.
It can therefore be useful to require a minimum assigned abundance before
considering a feature detected.

From **v1.1.0**, `--mincount` specifies the minimum number of assigned inserts
required for a feature to be retained. Features below the threshold are counted
as absent. The default value is `0`, which disables this filtering, and negative
values are not allowed.

With fractional multi-mapper sharing (`--multi=equal` or
`--multi=proportional`), assigned abundance and the amount removed by
`--mincount` can be fractional insert-equivalents.

The appropriate threshold depends on sequencing depth. For metagenomes or
metatranscriptomes with more than 10 million paired-end reads, we have
typically used a threshold of `10`.

## Physical inserts and profile insert-equivalents

The profiling summary distinguishes between **physical inserts** and
**insert-equivalents**.

The total, mapped, uniquely mapped, and multiply mapped counts describe
physical inserts in the alignment input. In contrast, the abundance profile
contains feature-level contributions. These contributions are reported as
insert-equivalents before sequence-length or unit normalization.

For `--multi=equal` and `--multi=proportional`, a successfully assigned
multi-mapping insert contributes fractions that sum to one insert-equivalent.
For `--multi=all`, however, every matched feature receives one full
insert-equivalent. The number of effective insert-equivalents can therefore
exceed the number of mapped or even total input inserts. This is expected
behavior for `--multi=all` and does not indicate that additional physical
inserts were observed.

When `--mincount` removes a feature, its feature-level contribution is counted
as a **Purged insert-equivalent**. The **Effective insert-equivalents** reported
in the header are calculated directly from the contributions remaining in the
retained profile.

When `--total` is supplied, contributions that are not represented by retained
features are included in the **Unknown** feature. This includes physically
unmapped inserts and, where applicable, contributions removed during
multi-mapper handling or by `--mincount`.

For example, suppose 1000 physical inserts were sequenced and 990 mapped.
Under `--multi=all`, multi-mapping could cause those 990 mapped inserts to
produce 1200 feature-level insert-equivalents. If `--mincount` subsequently
removes 50 insert-equivalents, the retained profile contains 1150 effective
insert-equivalents and:
```
Unknown = 10 unmapped inserts + 50 purged insert-equivalents = 60
```

The total pre-normalization profile mass is therefore:
```
1150 retained + 60 Unknown = 1210 insert-equivalents
```

This is the same profile mass as before `--mincount`:
```
1200 mapped feature contributions + 10 unmapped inserts = 1210
```

Thus, with `--multi=all`, profile mass is deliberately not constrained to the
number of physical inserts. `Unknown` should be interpreted as profile
contribution not represented by retained reference features, rather than
strictly as a count of physically unmapped inserts.

## Output format and useful information

Profile output is always gzip-compressed, irrespective of the filename used.

By default, the abundance table uses a pandas-compatible two-column header
containing `ID` and the sample label supplied with `--label`.
The `--pandas` option remains accepted for compatibility and requests the same
default format. Use `--no-pandas` to request the legacy output format, in which
the first column still contains feature IDs but its header label (`ID`) is omitted.

The header section of the output also includes commented provenance and
summary information about the profiling process, including the command used,
QNAME-grouping status, and insert-mapping statistics.

For example:

```text
# msamtools version: 1.2.0
# msamtools git commit: bedf810bcdaa
# Command line: msamtools profile --label test --unit rel --multi prop --total 110000 --genome genome.tsv -o test.profile.txt.gz alignments.bam
# QNAME grouping check: confirmed by input header SO:queryname
#
# Insert mapping statistics:
# Total inserts       : 110000 (100.00%)
# Mapped inserts      : 100000 ( 90.91%)
#   - Multiple mapped :   4423 (  4.02%)
#   - Uniquely mapped :  95577 ( 86.89%)
#
# Profile statistics:
# Purged insert-equivalents    :       0.00 (  0.00%) due to ambiguous mapping or low abundance features
# Effective insert-equivalents :  100000.00 ( 90.91%)
# Estimated sequence length for 'Unknown': 4228447 bp
#
ID      test
Unknown 0.09749158
Bacteroides_thetaiotaomicron__VPI-5482__GCF_000011065.1	0.18166258
Segatella_copri__HDE06__GCF_019249655.2	0.035903655
Mediterraneibacter_gnavus__ATCC_29149__GCF_009831375.1	0.045169526
Bacteroides_fragilis__NCTC_9343__GCF_000025985.1	0.072156601
Faecalibacterium_prausnitzii__APC918_95b__GCF_003312465.1	0.12621949
Roseburia_intestinalis__L1-82__GCF_900537995.1	0.06330019
Bifidobacterium_longum__NCC2705__GCF_000007525.1	0.090281838
Phocaeicola_vulgatus__ATCC_8482__GCF_000012825.1	0.1343378
Bifidobacterium_adolescentis__JCM_19861__GCF_036322895.1	0.027080085
Escherichia_coli__K-12_substr._MG1655__GCF_000005845.2	0.01805147
Escherichia_coli__Sakai_substr._RIMD_0509952__GCF_000008865.2	0.10834518
```

When sequence-length normalization is used and an **Unknown** fraction is
included, the average sequence length of the database is used as a proxy
sequence length for **Unknown**.

The sample label supplied through `--label` is included in the table header,
which makes it convenient to combine profiles from multiple samples.

## Combining multiple profiles

A standalone script for combining multiple msamtools profile files is available
as part of the [MIntO](https://github.com/arumugamlab/MIntO) pipeline:

[merge_profiles.R](https://github.com/arumugamlab/MIntO/blob/main/scripts/merge_profiles.R)

It can be used independently of the complete MIntO workflow and requires the
following R packages:

* `optparse`
* `R.utils`
* `data.table`
* `parallel`

## Command-line reference

The following reference is generated directly from:

```bash
./msamtools profile --help
```

<!-- BEGIN GENERATED HELP: profile -->
```text
Usage:
------

msamtools profile [-S] <bamfile> [--help] -o <file> --label=<string> [--genome=<string>] [--total=<int>] [--mincount=<int>] [--unit=<string>] [--pandas] [--no-pandas] [--nolen] [--multi=<string>]

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
  --total=<int>             number of high-quality inserts (mate-pairs/paired-ends) that were input to the aligner (default: unknown)
  --mincount=<int>          minimum number of inserts mapped to a feature, below which the feature is counted as absent (default: 0)
  --unit=<string>           unit of abundance to report {ab | rel | fpkm | tpm} (default: rel)
  --pandas                  use pandas-compatible two-column header: ID and sample label (default)
  --no-pandas               use legacy header with the first column unlabeled
  --nolen                   do not normalize the abundance (only relevant for ab or rel) for sequence length (default: normalize)
  --multi=<string>          how to deal with multi-mappers {all | equal | proportional | ignore} (default: proportional)

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

--total option:      In metagenomics, an unmapped insert could still be a valid
                     sequence, just missing in the database being mapped against.
                     This is the purpose of the '--total' option to track the
                     fraction of 'unknown' entities in the sample. Skipping --total
                     disables tracking the 'unknown' fraction. However, if the total
                     sequenced inserts were given, then there will be a new
                     feature added to denote the 'unknown' fraction.
Units of abundance:  Currently four different units are available.
                          'ab': raw insert-count abundance
                         'rel': relative abundance (default)
                        'fpkm': fragments per kilobase of sequence per million reads
                         'tpm': transcripts per million
                     If number of inserts input to the aligner is given via --total,
                     fpkm and tpm will behave differently than in RNAseq data,
                     as there is now a new entity called 'unknown'.
Alignment filtering: 'profile' expects that every alignment listed is considered
                     valid. For example, if one needs to filter alignments
                     based on alignment length, read length, alignment percent
                     identity, etc, this should have been done prior to
                     'profile'. Please see 'filter' for such filtering.
Multi-mapper inserts: Inserts mapping to multiple references need to be considered
                     carefully, as spurious mappings of promiscuous regions or
                     short homology could lead to incorrect abundances of
                     sequences. 'profile' offers four options for this purpose.
                     If an insert maps to N references at the same time:
                'ignore': insert is ignored.
                   'all': each reference gets 1 insert added.
                 'equal': each reference gets 1/N insert added.
          'proportional': each reference gets a fraction proportional to its
                          current assigned insert count; counts are initialized
                          from uniquely mapped inserts and refined iteratively.

```
<!-- END GENERATED HELP: profile -->
