# Compatibility

## Supported alignment formats

**msamtools** supports SAM and BAM input.

CRAM is currently unsupported and there are no plans to add CRAM support.

## CIGAR compatibility

Current develoption version of msamtools supports both traditional `M`
CIGAR operations and the extended `=` / `X` notation.

Where alignment matches and mismatches are relevant:

- `M` represents aligned query/reference positions whose match status must be
  determined from alignment metadata.
- `=` represents an explicit sequence match.
- `X` represents an explicit sequence mismatch.

Coverage treats `M`, `=` and `X` as aligned query bases contributing coverage
on reference positions. `D` and `N` advance along the reference but do not
contribute coverage.

## Compatibility with older msamtools versions

Releases up to and including v1.1.3 were primarily developed and tested with
alignments using traditional `M` CIGAR notation.

Some operations in those releases could produce incorrect results when input
alignments used extended `=` / `X` CIGAR notation. Alignments using
traditional `M` notation were not affected by this specific issue.

For reproducible analysis of alignments containing `=` or `X`, use msamtools
v1.2.0 or later.

### Impact on common metagenomic aligners

Common short-read aligners such as BWA-MEM, BWA-MEM2, and Bowtie2 use
traditional `M` CIGAR notation by default. Current Strobealign versions also
use `M` by default, although versions before 0.10.0 emitted `=`/`X` by
default. Bowtie2 and Strobealign can explicitly be configured to emit
extended `=`/`X` CIGAR notation.

## HTSlib and samtools

msamtools uses HTSlib internally for SAM/BAM I/O.

The HTSlib version used for a particular release is defined by the msamtools
build configuration. The accompanying msamtools + samtools container uses a
samtools executable matched to that HTSlib version.

Exact dependency and container versions for a release are given in the
corresponding release notes.
