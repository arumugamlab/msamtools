# Example workflows using msamtools

Similar to **samtools**, **msamtools** is designed to work on a stream,
avoiding creation of intermediate files. The examples below show common
workflows that combine msamtools with aligners and samtools.

> **NOTE:** Workflows using `--besthit` or `--uniqhit` require contiguous
> QNAME grouping. See [Input requirements](input-requirements.md#qname-grouping).

The human-read removal workflow is explained in more detail; the remaining
examples focus on the options most relevant to each task.

## Alignment and filtering in one step

If your aligner can write to `stdout`, its output can be piped directly to
**msamtools** and filtered on the fly.

### Task

Align **SAMPLE** (files `SAMPLE.1.fq.gz` and `SAMPLE.2.fq.gz`) to the
`bwa-mem2` database in `DB`; retain alignments with an alignment length of at
least `80 bp`, at least `95%` identity, and at least `80%` of the read aligned;
and write the output to `SAMPLE.DB.filtered.bam`.

### Command

```bash
bwa-mem2 mem DB SAMPLE.1.fq.gz SAMPLE.2.fq.gz \
   | msamtools filter -S -b -l 80 -p 95 -z 80 > SAMPLE.DB.filtered.bam
```

### Explanation

The command:

* aligns the reads with `bwa-mem2` and streams SAM output to `msamtools filter`;
* retains alignments with at least `80 bp` alignment length (`-l 80`), at least
  `95%` identity (`-p 95`), and at least `80%` of the read aligned (`-z 80`);
* writes the retained alignments in BAM format (`-b`).

## Removing human reads from human metagenomes

This workflow removes confidently mapped human reads while retaining unmapped
reads and short alignments that are more likely to represent spurious
host matches.

### Task

Align **SAMPLE** (files `SAMPLE.1.fq.gz` and `SAMPLE.2.fq.gz`) to the
human-genome `bwa-mem2` database in `HUMAN_DB`; trust alignments with at least
`30 bp` alignment length; and write the host-free reads as compressed FASTQ
files to `SAMPLE.hostfree.1.fq.gz` and `SAMPLE.hostfree.2.fq.gz`.

> **NOTE:** This is our standard workflow for host-sequence removal. Without
> the additional `30 bp` filtering step, more reads can be spuriously classified
> as host-derived. Passing the raw `bwa-mem2 mem` output directly to
> `samtools fastq` would therefore remove reads that this workflow retains.

### Command

```bash
bwa-mem2 mem HUMAN_DB SAMPLE.1.fq.gz SAMPLE.2.fq.gz \
  | msamtools filter -S -l 30 --invert --keep_unmapped -bu - \
  | samtools fastq -1 SAMPLE.hostfree.1.fq.gz -2 SAMPLE.hostfree.2.fq.gz -s /dev/null -0 /dev/null -c 6 -N -
```

### Explanation

The command:

* aligns the reads to the human reference with `bwa-mem2`, producing SAM output;
* pipes that SAM stream to `msamtools filter`;
* asks `msamtools filter` to:
    * read SAM input (`-S`);
    * define trusted host alignments as those with at least `30 bp` alignment
      length (`-l 30`);
    * invert that condition so alignments shorter than `30 bp` are retained
      (`--invert`);
    * retain unmapped reads as well (`--keep_unmapped`);
    * write an uncompressed BAM stream (`-bu`) for efficient piping;
* pipes the retained reads to `samtools fastq`;
* asks `samtools fastq` to:
    * write forward and reverse reads to separate files (`-1` and `-2`);
    * compress the FASTQ output (`-c 6`);
    * discard unpaired reads (`-s /dev/null -0 /dev/null`);
    * append `/1` and `/2` to read names (`-N`).

## Mapping a metagenome sample to a gene database and generating gene profiles

This workflow maps metagenomic reads to a gene catalog, keeps the best-scoring
alignments, and estimates relative gene abundances.

### Task

Align **SAMPLE** (files `SAMPLE.1.fq.gz` and `SAMPLE.2.fq.gz`) to the gene
catalog `bwa-mem2` database in `GENE_DB`; apply the same alignment thresholds
as in [Alignment and filtering in one step](#alignment-and-filtering-in-one-step);
retain only the highest-scoring hits; and write the gene profile to
`SAMPLE.profile.txt.gz`.

### Command

```bash
bwa-mem2 mem GENE_DB SAMPLE.1.fq.gz SAMPLE.2.fq.gz \
  | msamtools filter -S -bu -l 80 -p 95 -z 80 --besthit - \
  | msamtools profile --multi=proportional --label=SAMPLE --unit=rel -o SAMPLE.profile.txt.gz -
```

### Explanation

The command:

* filters alignments using the `80 bp`, `95%` identity, and `80%` aligned-query
  thresholds and retains the best-scoring hit(s) per read (`--besthit`);
* streams uncompressed BAM output directly to `msamtools profile`;
* shares multi-mapping inserts proportionally among reference genes
  (`--multi=proportional`);
* reports relative abundance (`--unit=rel`) using `SAMPLE` as the profile label.

## Estimating the number of inserts uniquely mapped to each gene in the database

This workflow is similar to the previous example, but retains only reads with a
unique best hit and reports insert counts rather than length-normalized relative
abundances.

### Task

Align **SAMPLE** (files `SAMPLE.1.fq.gz` and `SAMPLE.2.fq.gz`) to the gene
catalog `bwa-mem2` database in `GENE_DB`; apply the same alignment thresholds
as above; retain only reads with a unique highest-scoring hit; and write
insert counts for uniquely mapped inserts to `SAMPLE.profile.txt.gz`.

### Command

```bash
bwa-mem2 mem GENE_DB SAMPLE.1.fq.gz SAMPLE.2.fq.gz \
  | msamtools filter -S -bu -l 80 -p 95 -z 80 --uniqhit - \
  | msamtools profile --multi=ignore --label=SAMPLE --nolen --unit=ab -o SAMPLE.profile.txt.gz -
```

### Explanation

The command:

* retains only reads with a unique highest-scoring alignment (`--uniqhit`);
* ignores any remaining multi-mapping inserts during profiling
  (`--multi=ignore`);
* disables sequence-length normalization (`--nolen`);
* reports insert-count abundance (`--unit=ab`) using `SAMPLE` as the profile
  label.
