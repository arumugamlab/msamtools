# msamtools coverage

**coverage** program estimates per-position or fractional coverage of each sequence in
the BAM file.

## Per-position coverage of all sequences in BAM file <a name="pos-coverage"></a>

The per-position coverage output file is in the format of old Sanger quality files
with fasta headers and space-delimited numbers. As this can build up into
quite a large file, using the `-x` option will not print coverage for
sequences that did not have a single read mapped to them. Since their coverage
is essentially zero in each position, printing their coverage is just a
waste of space.

Here is an example per-position coverage command.
```bash
msamtools coverage -x -z -o sample1.coverage.txt.gz sample1.IGC.filtered.bam
```

## Fractional coverage of each sequence in BAM file <a name="frac-coverage"></a>

Sometimes it is useful to see which sequence from the BAM file has been observed in
the sample. And if yes, it is nice to know what fraction of the sequence has been
covered with alignments in the BAM file. For this one can use the `--summary`
option, which outputs fractional coverage and sequencing-coverage of each sequence.

Here is an example fractional coverage command.
```bash
msamtools coverage -z --summary -o sample1.coverage.summary.txt.gz sample1.IGC.filtered.bam
```

And here is an example output:
```text
cluster_001_consensus_length_3171293	0.05464883	0.25
cluster_002_consensus_length_2788722	0.99955930	10.79
cluster_003_consensus_length_6395848	0.99998921	38.10
cluster_004_consensus_length_2025181	0.99947906	31.14
cluster_005_consensus_length_3532514	0.99987346	70.04
```
First column names the sequence, 2nd column reports the fraction of that sequence that is covered
and the 3rd column gives sequencing-coverage. Apparently, the 5th genome has 70X coverage in
that sample!

A full description is given below:
```text
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
```
