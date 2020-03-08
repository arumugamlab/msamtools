# msamtools:  microbiome-related extension to samtools 

**msamtools** provides useful functions that are commonly used in microbiome
data analysis, especially when analyzing shotgun metagenomics 
or metatranscriptomics data.

# Table of Contents
 1. [Installation](#installation)   
  1.1. [System requirements](#sys-requirements)   
   1.1.1. [External packages](#required-packages)   
   1.1.2. [Required tools](#required-tools)   
  1.2. [For normal users](#normal-users)   
  1.3. [For advanced users](#advanced-users)   
   1.3.1. [Getting the source code](#source-code)   
    1.3.1.1. [Cloning the git repository](#git-clone)   
    1.3.1.2. [Downloading the ZIP file from github](#git-zip)   
    1.3.2. [Running autoconf and automake](#automake)   
    1.3.3. [Building the program](#build)   
 2. [msamtools](#msamtools)   
 3. [msamtools filter](#msamtools-filter)   
 4. [msamtools profile](#msamtools-profile)   
  4.1. [Keeping track of unmapped reads](#track-unmapped)   
 5. [msamtools coverage](#msamtools-coverage)   
 6. [msamtools summary](#msamtools-summary)   
## 1. Installation <a name="installation"></a>

### 1.1. System requirements <a name="sys-requirements"></a>

You should be able to use **msamtools** on any flavor of linux and UNIX. 
Although I have not tested it, it should also work on macOS. **msamtools** has
two dependencies (argtable2 and samtools), which it will automatically 
download and build. These two packages have their own requirements even though
I have tuned their configurations to a minimum. 

#### 1.1.1. External packages <a name="required-packages"></a>

The following packages are required to build **msamtools**:

 1. **zlib** development version
 2. **curses** development version (required by samtools)

Please make sure that these are installed in your system before trying to 
build.

#### 1.1.2. Required tools <a name="required-tools"></a>

While building **msamtools**, you will need some standard tools that are 
most likely installed in your system by default. I will still list them here
anyway to be sure:

 1. gcc
 2. gzip
 3. tar
 4. wget

If any of these is missing in your system, or cannot be found in your 
application path, please fix that first.

### 1.2. For normal users <a name="normal-users"></a>

If you are a normal user, then the easiest way is to obtain the package file
and build the program right away. The following commands were written when
version 0.9 was the latest, so please update the version number in the
commands below.

**Note:** Newer C compilers from gcc use `-std=gnu99` by default, which I had
not tested as my gcc version is quite outdated with `-std=gnu89` as default.
This leads to files not compiling when running `make` with new compilers. The
current fix for using the tarball for version 0.9 is to tell the compiler which
standard to use, using `CFLAGS="-std=gnu89"`. This extra option is a 
temporary fix only, and will not be needed in subsequent releases as the code 
gets upgraded to be compatible with `-std=gnu99`. 
(Thanks [Russel88](https://github.com/Russel88) for reporting this).

~~~
wget https://github.com/arumugamlab/msamtools/releases/download/0.9/msamtools-0.9.tar.gz
tar xfz msamtools-0.9.tar.gz
cd msamtools-0.9
./configure
make CFLAGS="-std=gnu89 -O2"
~~~

This should create `msamtools` executable amond others.

### 1.3. For advanced users <a name="advanced-users"></a>

If you are an advanced user who would like to contribute to the code base
or if you just like to do things the hard way, you can check out the source
code and build the program in a series of steps involving `autoconf` and
`automake`. If these names confuse you or scare you, then please follow the
instructions for [normal users](#normal-users).

#### 1.3.1. Getting the source code <a name="source-code"></a>

You can get **msamtools** code from github at 
<https://github.com/arumugamlab/msamtools>. 
You can either `git clone` it or download the ZIP file and extract the 
package.

##### 1.3.1.1. Cloning the git repository <a name="git-clone"></a>

You can get a clone of the repository if you wish to keep it up-to-date when
we release new versions or updates.

~~~
$ git clone https://github.com/arumugamlab/msamtools.git
Cloning into 'msamtools'...
remote: Enumerating objects: 41, done.
remote: Counting objects: 100% (41/41), done.
remote: Compressing objects: 100% (40/40), done.
remote: Total 41 (delta 0), reused 37 (delta 0), pack-reused 0
Unpacking objects: 100% (41/41), done.
$ cd msamtools-master
~~~

You can check the contents of the repository in *msamtools* directory.

##### 1.3.1.2. Downloading the ZIP file from github <a name="git-zip"></a>

You can download the repository's snapshot as on the day of download by:
~~~
$ wget https://github.com/arumugamlab/msamtools/archive/master.zip
--2020-02-10 16:45:39--  https://github.com/arumugamlab/msamtools/archive/master.zip
Resolving github.com (github.com)... 140.82.118.4
Connecting to github.com (github.com)|140.82.118.4|:443... connected.
HTTP request sent, awaiting response... 302 Found
Location: https://codeload.github.com/arumugamlab/msamtools/zip/master [following]
--2020-02-10 16:45:39--  https://codeload.github.com/arumugamlab/msamtools/zip/master
Resolving codeload.github.com (codeload.github.com)... 192.30.253.120
Connecting to codeload.github.com (codeload.github.com)|192.30.253.120|:443... connected.
HTTP request sent, awaiting response... 200 OK
Length: unspecified [application/zip]
Saving to: 'master.zip'

    [ <=>                                   ] 55,635       280K/s   in 0.2s
$ unzip master.zip
$ cd msamtools-master
~~~

##### 1.3.2. Running autoconf and automake <a name="automake"></a>

You can check the contents of the repository in the package directory.
~~~
$ ls
configure.ac    mBamVector.h    mMatrix.h        msamtools.c  splitSeq.c
deps            mCommon.c       msam_coverage.c  mSequence.c  versions.txt
filterSeq.c     mCommon.h       msam_filter.c    mSequence.h  zoeTools.c
LICENSE         mCompress.c     msam.h           mSimulate.c  zoeTools.h
Makefile.am     mCompress.h     msam_helper.c    mSimulate.h
make_readme.sh  mDefinitions.h  msam_profile.c   README.md
mBamVector.c    mMatrix.c       msam_summary.c   seqUtils.c
~~~

You will note that the `configure` script does not exist in the package.
This is because you need to generate the `configure` script using 
`aclocal`, `autoconf` and `automake`.
~~~
aclocal
autoconf
mkdir build-aux
automake --add-missing
~~~

##### 1.3.3. Building the program <a name="build"></a>

You can then build msamtools as follows:
~~~
./configure
make
~~~

## 2. msamtools <a name="msamtools"></a>

This is the master program that you call with the subprogram options. There 
are currently 4 subprograms that you can call as shown below.
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

These represent the different analysis of SAM/BAM files you can perform using
**msamtools**.

## 3. msamtools filter <a name="msamtools-filter"></a>

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

## 4. msamtools profile <a name="msamtools-profile"></a>

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
msamtools profile --multi=proportional --label=sample1 -z -o sample1.IGC.profile.txt.gz sample1.IGC.filtered.bam
~~~
The above command estimates relative abundance of IGC genes after sharing
multi-mapper reads proportionately between the genes (see below).

In the spirit of **samtools** programs, **msamtools** programs can also 
stream between each other. Therefore, a single command to **filter** and **profile** 
would look like:
~~~
msamtools filter -b -u -l 80 -p 95 -z 80 --besthit sample1.IGC.bam | msamtools profile --multi=proportional --label=sample1 -o sample1.IGC.profile.txt.gz -z -
~~~

### 4.1. Keeping track of unmapped reads <a name="track-unmapped"></a>

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

If you would like to include the Unknown fraction in the profile, you should
tell the profiler how many reads were originally sequenced. Then the profiler will estimate
the unmapped reads based on how many reads are in the bam/sam file, and then use
it in the profiling stage.

~~~
 # Get number of entries in the fwd fastq file = number of inserts
lines=$(zcat sample1.1.fq.gz | wc -l)
entries=$(expr $lines / 4)   # There are 4 lines per entry

 # Use total reads in profiler
msamtools filter -b -u -l 80 -p 95 -z 80 --besthit sample1.IGC.bam | msamtools profile --multi=proportional --label=sample1 -o sample1.IGC.profile.txt.gz -z --total=$entries -
~~~

A full description is given below:
~~~
Usage:
------

msamtools profile [-Sz] <bamfile> [--help] -o <file> --label=<string> [--multi=<string>] [--total=<int>]

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
  --total=<int>             number of high-quality inserts (mate-pairs/paired-ends) that were input to the aligner (default: 0)
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

## 5. msamtools coverage <a name="msamtools-coverage"></a>

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

## 6. msamtools summary <a name="msamtools-summary"></a>

**summary** program summarizes alignments given in the BAM file. It can
also provide distributions of certain features across alignments.

Here is an example summary command that reports all alignments:
~~~
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

~~~

Here is another example summary command that reports the distribution
of unmapped bases:
~~~
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
~~~
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

