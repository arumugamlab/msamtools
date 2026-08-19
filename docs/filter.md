# msamtools filter

The **filter** command provides alignment filtering based on percent identity,
alignment length, aligned fraction of read length, or combinations thereof.
For example, when mapping metagenomic reads to a database for species-level
annotation, we typically discard alignments below 95% sequence identity.

Here is an example filtering command one would use after mapping metagenomic
reads to the Integrated Gene Catalog (IGC) consisting of 9.9 million genes
(Li *et al.*, **Nat. Biotech.** 2014).

```bash
msamtools filter -b -l 80 -p 95 -z 80 --besthit sample1.IGC.bam > sample1.IGC.filtered.bam
```

The command selects alignments that are at least 80 bp long, have at least
95% sequence identity, and cover at least 80% of the read. Additionally,
`--besthit` retains the highest-scoring alignment(s) for each read. Multiple
alignments with the same highest alignment score are retained.

A full description is given below:

```text
Usage:
-------

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
  -p <int>                  min. sequence identity of alignment, in percentage, integer between 0 and 100; requires MD or NM field to be present (default: 0)
  --ppt=<int>               min/max sequence identity of alignment, in parts per thousand, integer between -1000 and 1000; requires MD or NM field to be present (default: 0)
                            NOTE:
                            -----
                                  When using --ppt, +ve values mean minimum ppt and -ve values mean maximum ppt.
                                  E.g., '--ppt 950' will report alignments with ppt >= 950,
                                  and '--ppt -950' will report alignments with ppt <= 950.
  -z <int>                  min. percent of the query that must be aligned, between 0 and 100 (default: 0)
  -k, --keep_unmapped       report unmapped reads, when filtering using upper-limit thresholds (default: false)
  -v, --invert              invert the effect of the filter (default: false)
                            CAUTION:
                            --------
                                  When using --invert or -v, be precise in what needs to be inverted.
                                  Adding -v gives you the complement of what you get without -v.
                                  Sometimes, this might be counter-intuitive.
                                  E.g., '-l 65 -p 95' will report alignments that are (>=65bp AND >=95%).
                                        '-l 65 -p 95 -v' will not report (<65bp AND <95%), as one might think.
                                        '-l 65 -p 95 -v' will report NOT(>=65bp AND >=95%) which is (<65bp OR <95%).
                                        Notice the 'OR' in the final logical operation. This means that
                                        an alignment that fails one condition will still be reported if it
                                        satisfies the other condition.
                                        If you only wanted alignments that are below 95%, then specify '-p 95 -v'
  --rescore                 rescore alignments using MD or NM fields, in that order (default: false)

Special filters:
----------------

The following special filters cannot be combined with -v, but require:
  (1) the alignments to be sorted by name,
  (2) AS field (alignment score) to be present, unless --rescore is used.
You can usually achieve sorting by:
  samtools sort -n -T tmp.sort input.bam | msamtools filter --besthit -
If AS is missing, you can rescore alignments by:
  samtools sort -n -T tmp.sort input.bam | msamtools filter --rescore --besthit -

  --besthit                 keep all highest scoring hit(s) per read (default: false)
  --uniqhit                 keep only one highest scoring hit per read, only if it is unique (default: false)
```

The `filter` command operates at the read level. With paired-end data,
`--besthit` and `--uniqhit` select alignments independently for READ1 and
READ2. Output for each QNAME group is normalized so that retained READ1
alignments are written before retained READ2 alignments.

When using `--besthit` or `--uniqhit`, all alignments with the same **QNAME**
must occur as one contiguous group in the input. Alignments belonging to
READ1 and READ2 may be interleaved within that group; the mates do not need
to occur as separate blocks.

Typical read-mapper output may already satisfy this requirement. If the BAM
file has subsequently been processed or reordered, sort it by QNAME before
filtering, for example:

```bash
samtools sort -n input.bam | msamtools filter --besthit -
```
