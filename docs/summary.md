# msamtools summary

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
```
