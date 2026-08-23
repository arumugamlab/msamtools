# Input requirements

**msamtools** reads SAM or BAM alignment files. SAM input should be indicated
with `-S`; BAM input is used by default.

Some msamtools operations process each alignment independently and do not require
any particular record ordering. Others need all alignments belonging to the
same read or insert to occur together.

## QNAME grouping

For operations that work across multiple alignments of the same read or
insert, all records sharing a **QNAME** must occur as one contiguous group in
the input.

A full lexical sort by QNAME is not inherently required; the essential
requirement is that a QNAME must not reappear after records from another QNAME
have begun.

QNAME-sorted input is therefore the usual and recommended way to satisfy this
requirement:

```bash
samtools sort -n input.bam -o input.qname.bam
```

When streaming:

```bash
samtools sort -n input.bam | msamtools <command> ... -
```

For paired-end data, READ1 and READ2 records belonging to the same QNAME may
be interleaved within the QNAME group. They do not need to occur as separate
mate-specific blocks.

## Commands requiring contiguous QNAME groups

Contiguous QNAME grouping is required for:

- `filter --besthit`
- `filter --uniqhit`
- `profile`
- `summary --count`

For `filter --besthit` and `filter --uniqhit`, alignments are selected
independently for READ1 and READ2 within each QNAME group.

The `profile` command uses QNAME groups to identify inserts and correctly
handle inserts mapping to multiple references.

The `summary --count` option counts mapped inserts represented by contiguous
QNAME groups.

## Coordinate-sorted BAM files

Coordinate-sorted BAM files generally do not satisfy the QNAME-grouping
requirement because alignments from the same read can occur at different
positions in the file.

Before using an operation that requires QNAME grouping, name-sort such files,
for example:

```bash
samtools sort -n input.coord-sorted.bam \
    | msamtools profile --label sample -o sample.profile.txt.gz -
```

Operations that process each alignment independently do not require
QNAME-sorted input.

## SAM and BAM input

BAM is the default input format.

For SAM input, use `-S` where supported:

```bash
msamtools filter -S input.sam
```

msamtools commands that accept alignment input also accept `-` as input,
allowing SAM/BAM data to be streamed from another program.

## Alignment metadata requirements

Some filtering operations depend on optional SAM alignment tags.

Percent-identity filtering with `-p` or `--ppt` requires alignment mismatch
information from either the `MD` or `NM` field.

`filter --besthit` and `filter --uniqhit` normally use the `AS` alignment-score
field. If `AS` is unavailable, `--rescore` can be used to recompute alignment
scores from `MD` or `NM` information.

For command-specific requirements and options, see the documentation for
[`filter`](filter.md), [`profile`](profile.md), and [`summary`](summary.md).
