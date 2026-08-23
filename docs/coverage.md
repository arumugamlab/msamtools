# msamtools coverage

The **coverage** command estimates per-position or per-sequence coverage of each
reference sequence in the SAM/BAM file.

Coverage is defined as **aligned query-base depth on reference positions**.
Thus CIGAR operations `M`, `=` and `X` contribute coverage, while `D` and `N`
advance along the reference but do not contribute coverage.

Output is always gzip-compressed. The filename is used exactly as provided;
msamtools does not automatically append `.gz`.

## Per-position coverage of all sequences

Per-position coverage is written in a format similar to old Sanger quality
files, with FASTA-style headers followed by space-delimited coverage values.

For large reference databases in which only a small fraction of reference
sequences are covered, `-x` (`--skipuncovered`) avoids writing long arrays of
zeros for completely uncovered sequences.

For example:

```bash
msamtools coverage -x -o sample1.coverage.txt.gz sample1.IGC.filtered.bam
```

By default, uncovered reference sequences are also reported, with zero
coverage at every position.

## Coverage summary for each sequence

The `--summary` option reports one line per reference sequence rather than
per-position coverage.

For example:

```bash
msamtools coverage --summary -o sample1.coverage.summary.txt.gz sample1.IGC.filtered.bam
```

Example output:

```text
cluster_001_consensus_length_3171293	0.05464883	0.25
cluster_002_consensus_length_2788722	0.99955930	10.79
cluster_003_consensus_length_6395848	0.99998921	38.10
cluster_004_consensus_length_2025181	0.99947906	31.14
cluster_005_consensus_length_3532514	0.99987346	70.04
```

The columns are:

1. reference sequence name
2. fraction of reference positions covered by at least one aligned query base
3. mean aligned query-base depth across the complete reference sequence

Thus, in the example above, the fifth reference sequence has approximately
70× mean coverage.

The `-x` (`--skipuncovered`) option can also be used with `--summary` to omit
reference sequences with no covered positions.

## Command-line reference

The following reference is generated directly from:

```bash
./msamtools coverage --help
```

<!-- BEGIN GENERATED HELP: coverage -->
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
  -z, --gzip                compress output file using gzip (default; option retained for backward compatibility)

Description:
------------

Produces per-position sequence coverage information for all reference sequences
in the BAM file. Output is similar to old-style quality files from the Sanger
sequencing era, with a fasta-style header followed by lines of space-delimited
numbers.

Coverage is defined as aligned query-base depth on reference positions. Thus
CIGAR operations 'M', '=' and 'X' contribute coverage, while 'D' and 'N' do not.

For large reference databases, option '-x' comes in handy when only a small
fraction of reference sequences are covered.

Output is always gzip-compressed, but file names do NOT automatically get a
'.gz' extension; specify the desired full output file name.
```
<!-- END GENERATED HELP: coverage -->
