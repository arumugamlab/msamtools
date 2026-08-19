# msamtools profile

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
```bash
msamtools profile --multi=proportional --label=sample1 --unit=rel -o sample1.IGC.profile.txt.gz sample1.IGC.filtered.bam
```
The above command estimates relative abundance of IGC genes after sharing
multi-mapper reads proportionately between the genes (see below).

In the spirit of **samtools** programs, **msamtools** programs can also
stream between each other. Therefore, a single command to **filter** and **profile**
would look like:
```bash
msamtools filter -b -u -l 80 -p 95 -z 80 --besthit sample1.IGC.bam \
  | msamtools profile --multi=proportional --label=sample1 --unit=rel -o sample1.IGC.profile.txt.gz -
```

or for mapping to scaffolds that are grouped in metagenome-assembled MAGs using:
```bash
msamtools filter -b -u -l 80 -p 95 -z 80 --besthit sample1.myMAGs.bam \
  | msamtools profile --multi=proportional --label=sample1 --unit=rel --genome myMAGs.genome.def -o sample1.myMAGs.profile.txt.gz -
```

## Profiling genomes or MAGs <a name="profiling-genomes"></a>

Starting from **v1.0.0**, **profile** program supports profiling of genomes defined by a set of
sequences. This requires a tab-delimited definition file of the following format:
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

If this information is stored in a file called `myMAGs.genome.def`, then you can
run the profiler as follows to get profiles at the genome level.
```bash
msamtools filter -b -u -l 80 -p 95 -z 80 --besthit sample1.myMAGs.bam \
  | msamtools profile --multi=proportional --label=sample1 --unit=rel --genome myMAGs.genome.def -o sample1.myMAGs.profile.txt.gz -
```

## Units of abundance <a name="abundance-units"></a>

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

## Keeping track of unmapped reads <a name="track-unmapped"></a>

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

```bash
# Get number of entries in the fwd fastq file = number of inserts
lines=$(zcat sample1.1.fq.gz | wc -l)
entries=$(expr \$lines / 4)   # There are 4 lines per fastq entry

# Use total reads in profiler
msamtools filter -b -u -l 80 -p 95 -z 80 --besthit sample1.IGC.bam \
  | msamtools profile --multi=proportional --label=sample1 -o sample1.IGC.profile.txt.gz --total=$entries -
```

## Avoiding extremely-low-abundant features <a name="profile-mincount"></a>

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

## Useful information in the output file <a name="profile-output"></a>

The header section of the output file includes a few lines of comment that
are hopefully useful. Here is an example:

```text
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
```

The commented lines are self-explanatory, and could be useful in getting
quick summary of the profiling process. Since length-normalization is not
turned off, the average sequence length for the entire database is used
as a proxy sequence length for the **Unknown** fraction.

The first line includes the name of the sample provided via `--label`.
This is for conveniently combining output from multiple files. A script that
combines output does not need external information to create a table with the
right sample name in a row/column.

A full description is given below:
```text
Usage:
------

msamtools profile [-S] <bamfile> [--help] -o <file> --label=<string> [--genome=<string>] [--total=<int>] [--mincount=<int>] [--unit=<string>] [--pandas] [--nolen] [--multi=<string>]

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
  --pandas                  print two columns (ID, sample-label) as header compatible with python pandas (default: only sample label)
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

--total option:      In metagenomics, an unmapped insert could still be a valid
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
                     sequences. 'profile' offers three options for this purpose.
                     If an insert maps to N references at the same time:
                'ignore': insert is ignored.
                   'all': each reference gets 1 insert added.
                 'equal': each reference gets 1/N insert added.
          'proportional': each reference gets a fraction proportional to its 
                          reference-sequence-length-normalized relative 
                          abundance estimated only based on uniquely
                          mapped inserts.
```
