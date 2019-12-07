# msamtools:  microbiome-related extension to samtools

**msamtools** provides useful functions that are commonly used in microbiome
data analysis, especially when analyzing shotgun metagenomics 
or metatranscriptomics data.

## msamtools

This is the master program that you call with the subprogram options.
~~~

msamtools: Metagenomics-related extension to samtools.
  Version: 0.9

Usage:
------

msamtools <command> [options]

Command:  filter         filter alignments based on alignment statistics
          profile        estimate relative abundance profile of reference sequences in bam file
          coverage       estimate per-base read coverage of each reference sequence
          summary        summarize alignment statistics per read in a table format
~~~

## msamtools filter

**filter** program provides alignment filtering based on percent identity,
read length, aligned fraction of read length, or combinations thereof.
For example, in mapping metagenomic reads to a database for species-level 
annotation, we typically throw out alignments <95% sequence identity.

Here is an example filtering command one would use after mapping metagenomic
reads to the Integrated Gene Catalog consisting 9.9 million genes 
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

msamtools filter [-buhSv] <bamfile> [--help] [-l <int>] [-p <int>] [--ppt=<int>] [-z <int>] [--rescore] [--besthit] [--uniqhit]

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

## msamtools profile

**profile** program provides sequence abundance profiling functionality.
Relative abundance of each sequence in the BAM file is reported. Abundance
of a sequence is estimated as number of inserts (mate-pairs/paired-ends) 
mapped to that sequence, divided by its length. Reads mapping to multiple 
sequences can be shared across the sequences in three different ways 
(please see below).  Finally, relative abundance
is estimated by normalizing the total abundances across the BAM file.

We highly recommend that you filter the alignments before sending to the
**profile** program, as it considers each alignment to be important (it 
does not look at alignment quality, for example).

Here is an example profiling command one would use after mapping metagenomic
reads to IGC.
~~~
msamtools profile --multi=proportional --label=sample1 -z -o sample1.profile.txt.gz sample1.IGC.filtered.bam
~~~
The above command estimates relative abundance of IGC genes after sharing
multi-mapper reads proportionately between the genes (see below).

In the spirit of **samtools** programs, **msamtools** programs can also 
stream between each other. Therefore, a single command to **filter** and **profile** 
would look like:
~~~
msamtools filter -b -u -l 80 -p 95 -z 80 --besthit sample1.IGC.bam | msamtools profile --multi=proportional --label=sample1 -o sample1.profile.txt.gz -z -
~~~

#### Keeping track of unmapped reads

By default, **profile** command will generate relative abundances that sum to `100%` or `1` 
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

If you do not have the number of unmapped reads handy, here is one way to
estimate it using **summary** command (see Section 
**msamtools summary** for more details):

~~~
# Get number of entries in the fwd fastq file = number of inserts
lines=$(zcat sample1.1.fq.gz | wc -l)
entries=$(expr $lines / 4)   # There are 4 lines per entry

# Get number of inserts in the BAM file
mapped=$(msamtools summary --count sample1.IGC.filtered.bam)

# Get unmapped
unmapped=$(expr $entries - $mapped)

# Use unmapped in profiler
msamtools profile --label sample1 -o sample1.IGC.profile --unmapped=$unmapped sample1.IGC.filtered.bam
~~~

A full description is given below:
~~~
Usage:
------

msamtools profile [-Sz] <bamfile> [--help] -o <file> --label=<string> [--multi=<string>] [--unmapped=<int>]

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
  --multi=<string>          how to deal with multi-mappers {all | equal | proportional} (default: proportional)
  --unmapped=<int>          number of inserts (mate-pairs/paired-ends) that were unmapped (default: 0)
  -z, --gzip                compress output file using gzip (default: false)

Description
-----------

Produces a relative abundance profile of all reference sequences in a BAM file
based on the number of read-pairs (inserts) mapping to each reference sequence.
It can work with genome-scale reference sequences while mapping to a database 
of sequenced genomes, but can also work with gene-scale sequences such as in the
Integrated Gene Catalog from human gut microbiome (Li et al, Nat biotech 2014).

In the output file, each sequence in the BAM file gets a line with its abundance.
They are presented in the order in which they appear in the BAM header. <label>
is used as the first line, so that reading or 'joining' these files is easier.

If using '-z', output file does NOT automatically get '.gz' extension. This is 
up to the user to specify the correct full output file name.

Alignment filtering: It expects that every alignment listed is considered 
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

## msamtools coverage

**coverage** program estimates per-position coverage of each sequence in
the BAM file. The output file is in the format of old Sanger quality files
with fasta headers and space-delimited numbers. As this can build up into
quite a large file, using the `-x` option will not print coverage for
sequences that did not have a single read mapped to them. Since their coverage
is essentially zero in each position, printing their coverage is just a 
waste of space.

Here is an example coverage command.
~~~
msamtools coverage -x -z -o sample1.coverage.txt.gz sample1.IGC.filtered.bam
~~~

A full description is given below:
~~~
Usage:
------

msamtools coverage [-Sxz] <bamfile> [--help] -o <file> [-w <int>]

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
  -x, --skipuncovered       do not report coverage for sequences without aligned reads (default: false)
  -w, --wordsize=<int>      number of words (coverage values) per line (default: 17)
  -z, --gzip                compress output file using gzip (default: false)

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

## msamtools summary

**summary** program summarizes alignments given in the BAM file. It can
also provide distributions of certain features across alignments.

Here is an example summary command that reports all alignments:
~~~
msamtools summary sample1.IGC.bam
~~~

Here is another example summary command that reports the distribution
of mapped bases:
~~~
msamtools summary --stats=mapped sample1.IGC.bam
~~~

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

