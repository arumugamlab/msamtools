# msamtools summary

The **summary** command reports alignment-level statistics from a SAM/BAM file.
It can also report distributions of selected alignment statistics or count
mapped inserts/QNAME groups.

By default, `summary` reports one tab-delimited line for each mapped,
non-secondary alignment. Supplementary alignments are retained.

Here is an example:

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
ERR688505.5_FCD1R4CACXX:7:1101:1765:2211#CAGCGGCG       78      MH0188_GL0130879        78      78      100.0
```

The six fields are:

1. QNAME
2. query length
3. reference sequence name
4. glocal alignment length
5. number of matches
6. percent identity

The glocal alignment length includes unaligned query bases, corresponding to a
global alignment with respect to the query and a local alignment with respect
to the reference.

## Alignment-statistic distributions

The `--stats` option reports a distribution rather than one line per alignment.

Available statistics are:

- `mapped` — number of matched query bases
- `unmapped` — number of query bases that are not matches
- `edit` — edit distance
- `score` — matches minus edit distance

As with the default summary output, unmapped and secondary alignments are
excluded.

For example:

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
...
```

Here, 50,033 alignments had zero query bases that were not matches, 24,425 had
one, 12,495 had two, and so on.

## Excluding alignments near reference edges

The `-e` / `--edge` option excludes alignments that occur too close to either
end of the reference sequence.

For example:

```bash
msamtools summary -e 10 input.bam
```

retains an alignment only when at least 10 reference bases remain between the
alignment and each end of the reference sequence.

`-e` must be non-negative.

## Counting mapped inserts

The `-c` / `--count` option reports the number of mapped inserts, represented
as distinct contiguous QNAME groups in the input.

```bash
msamtools summary --count input.bam
```

Unmapped records are not counted.

Because counting is based on QNAME groups, records belonging to the same QNAME
must occur contiguously in the input. QNAME-sorted input is therefore the
appropriate input organization when using `--count`.

`--count` cannot be combined with either `--stats` or `-e`.

## Command-line reference

```text
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
  -c, --count               count number of mapped inserts/QNAME groups in BAM file (default: false)
  --stats=<string>          {mapped|unmapped|edit|score} only report readcount distribution for specified stats, not read-level stats (default: none)

Description
-----------
Prints summary of alignments in the given BAM/SAM file. By default, it prints
a summary line per alignment entry in the file. The summary is a tab-delimited
line with the following fields:
        qname, query_length, target_name, glocal_align_len, matches, percent_identity
glocal_align_len includes the unaligned qlen mimicing a global alignment
in the query and local alignment in target, thus glocal.

With --stats option, summary is consolidated as distribution of read counts
for a given measure.
   --stats=mapped   - distribution for number of mapped query bases
   --stats=unmapped - distribution for number of unmapped query bases
   --stats=edit     - distribution for edit distances
   --stats=score    - distribution for score=match-edit
```
