#!/bin/bash

cat <<EOF
# msamtools:  microbiome-related extension to samtools

**msamtools** provides useful functions that are commonly used in microbiome
data analysis, especially when analyzing shotgun metagenomics 
or metatranscriptomics data.

## msamtools

This is the master program that you call with the subprogram options.
~~~
`./msamtools help`
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
`./msamtools filter --help`
~~~

It is important to note that **filter** command works at the read level.
E.g., using \`--besthit\` will independently choose the best-scoring alignments
for the forward and reverse reads.

Additionally, when using \`--besthit\` or \`--uniqhit\`, the input BAM file
is expected to be sorted by **QNAME**. Each mate should occur as a continuous
group followed by the other mate. Typical output from read-mappers do
produce output that follows this. However, if you had processed the BAM 
files yourself and modified things, it could affect output from **filter**.
In such cases, you should sort them again using \`samtools sort -n\` so they
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

By default, **profile** command will generate relative abundances that sum to \`100%\` or \`1\` 
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
lines=\$(zcat sample1.1.fq.gz | wc -l)
entries=\$(expr \$lines / 4)   # There are 4 lines per entry

# Get number of inserts in the BAM file
mapped=\$(msamtools summary --count sample1.IGC.filtered.bam)

# Get unmapped
unmapped=\$(expr \$entries - \$mapped)

# Use unmapped in profiler
msamtools profile --label sample1 -o sample1.IGC.profile --unmapped=\$unmapped sample1.IGC.filtered.bam
~~~

A full description is given below:
~~~
`./msamtools profile --help`
~~~

## msamtools coverage

**coverage** program estimates per-position coverage of each sequence in
the BAM file. The output file is in the format of old Sanger quality files
with fasta headers and space-delimited numbers. As this can build up into
quite a large file, using the \`-x\` option will not print coverage for
sequences that did not have a single read mapped to them. Since their coverage
is essentially zero in each position, printing their coverage is just a 
waste of space.

Here is an example coverage command.
~~~
msamtools coverage -x -z -o sample1.coverage.txt.gz sample1.IGC.filtered.bam
~~~

A full description is given below:
~~~
`./msamtools coverage --help`
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
`./msamtools summary --help`
~~~

EOF
