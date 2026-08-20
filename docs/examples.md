# Example workflows using msamtools

Similar to **samtools**, **msamtools** is designed to work on a stream,
avoiding creation of intermediate files. Here are some example
workflows using streams that **msamtools** will be useful in.

> **NOTE:** Workflows using `--besthit` or `--uniqhit` require contiguous
> QNAME grouping. See [Input requirements](input-requirements.md#qname-grouping).

## Alignment and filtering in one step

If your aligner can write to `stdout`, then you can directly pipe the output
to **msamtools** and filter on the fly.

### Task
Align **SAMPLE** (files `SAMPLE.1.fq.gz` and `SAMPLE.2.fq.gz`) to the
`bwa-mem2` database in `DB`; retain alignments at least `80bp` long,
with at least `95%` identity, covering at least `80%` of the read;
and write output to `SAMPLE.DB.filtered.bam`.

### Command
```bash
bwa-mem2 mem DB SAMPLE.1.fq.gz SAMPLE.2.fq.gz \
   | msamtools filter -S -b -l 80 -p 95 -z 80 > SAMPLE.DB.filtered.bam
```

### Explanation

The command above

* aligns using `bwa-mem2` that generates `SAM` format
* pipes the output to `msamtools`
* asks `msamtools filter` to
    * read `SAM` format (`-S`)
    * filter alignments that are
        * at least `80bp` long (`-l 80`)
        * at least `95%` identity (`-p 95`)
        * at least `80%` of the read aligned (`-z 80`)
    * write output in `BAM` format (`-b`)

## Removing human reads from human metagenomes

Here is an example workflow to filter human reads using `bwa-mem2`.

### Task
Align **SAMPLE** (files `SAMPLE.1.fq.gz` and `SAMPLE.2.fq.gz`) to the
human genome `bwa-mem2` database in `HUMAN_DB`; trust alignments with
at least `30bp` alignment length; and write the host-free reads as compressed
`fastq` files to `SAMPLE.hostfree.1.fq.gz` and `SAMPLE.hostfree.2.fq.gz`.

>**Note:** This is our standard workflow for host-sequence removal.
Without the extra filtering by `30bp`, we lose significantly more
reads that are spuriously flagged as host-derived. Directly going
from `bwa-mem2 mem` output to `samtools` is not advisable, as this
will remove useful reads from your sample.

### Command
```bash
bwa-mem2 mem HUMAN_DB SAMPLE.1.fq.gz SAMPLE.2.fq.gz \
  | msamtools filter -S -l 30 --invert --keep_unmapped -bu - \
  | samtools fastq -1 SAMPLE.hostfree.1.fq.gz -2 SAMPLE.hostfree.2.fq.gz -s /dev/null -0 /dev/null -c 6 -N -
```

### Explanation

The command above

* aligns using `bwa-mem2` that generates `SAM` format
* pipes the output to `msamtools`
* asks `msamtools filter` to get reads that are not human by
    * reading `SAM` format (`-S`)
    * filtering alignments that are at least `30bp` long (`-l 30`)
    * negating that condition and retaining alignments below `30bp` (`--invert`)
    * while retaining also the unmapped reads (`--keep_unmapped`)
    * writing output in uncompressed `BAM` format (`-bu`)
* then pipes the output to `samtools`
* asks `samtools` to make `fastq` files
    * write compressed forms (`-c 6`)
    * of forward and reverse reads to separate files (`-1` and `-2`)
    * while ignoring unpaired reads (`-s /dev/null -0 /dev/null`)
    * and appending `/1` and `/2` to the reads (`-N`)

## Mapping a metagenome sample to a gene database and generating gene profiles

Here is an example workflow one would use after mapping metagenomic reads to IGC.

### Task

Align **SAMPLE** (files `SAMPLE.1.fq.gz` and `SAMPLE.2.fq.gz`) to the
gene catalog `bwa-mem2` database in `GENE_DB`; filter as in
[Alignment and filtering in one step](#alignment-and-filtering-in-one-step)
but retain only the highest-scoring hits; and write profile of all
genes to `SAMPLE.profile.txt.gz`.

### Command
```bash
bwa-mem2 mem GENE_DB SAMPLE.1.fq.gz SAMPLE.2.fq.gz \
  | msamtools filter -S -bu -l 80 -p 95 -z 80 --besthit - \
  | msamtools profile --multi=proportional --label=SAMPLE --unit=rel -o SAMPLE.profile.txt.gz -
```

### Explanation

The command above

* aligns using `bwa-mem2` that generates `SAM` format
* then pipes the output to `msamtools filter` to
    * read `SAM` format (`-S`)
    * filter alignments that are
        * at least `80bp` long (`-l 80`)
        * at least `95%` identity (`-p 95`)
        * at least `80%` of the read aligned (`-z 80`)
    * keep only best-scoring hits per read (`--besthit`)
    * write output in uncompressed `BAM` format (`-bu`)
* then pipes the output to `msamtools profile` to
    * share multihit inserts proportionally among hits (`--multi=proportional`)
    * calculate relative abundance profiles (`--unit=rel`)
    * use sample label `SAMPLE` in the output file (`--label=SAMPLE`)
    * and write compressed output to `SAMPLE.profile.txt.gz` (`-o`)

## Estimating the number of inserts uniquely mapped to each gene in database

Here is an example workflow one would use after mapping metagenomic reads to IGC.

### Task

Align **SAMPLE** (fastq files `SAMPLE.1.fq.gz` and `SAMPLE.2.fq.gz`) to the gene catalog `bwa-mem2` database in `GENE_DB`; filter as in [Mapping a metagenome sample to a gene database and generating gene profiles](#mapping-a-metagenome-sample-to-a-gene-database-and-generating-gene-profiles) but retain only reads uniquely mapping to a single reference; and write insert-counts for uniquely mapped inserts for all genes to `SAMPLE.profile.txt.gz`.

### Command
```bash
bwa-mem2 mem GENE_DB SAMPLE.1.fq.gz SAMPLE.2.fq.gz \
  | msamtools filter -S -bu -l 80 -p 95 -z 80 --uniqhit - \
  | msamtools profile --multi=ignore --label=SAMPLE --nolen --unit=ab -o SAMPLE.profile.txt.gz -
```

### Explanation

The command above

* aligns using `bwa-mem2` that generates `SAM` format
* then pipes the output to `msamtools filter` to
    * read `SAM` format (`-S`)
    * filter alignments that are
        * at least `80bp` long (`-l 80`)
        * at least `95%` identity (`-p 95`)
        * at least `80%` of the read aligned (`-z 80`)
    * keep only best-scoring hits per read provided best-hit is unique (`--uniqhit`)
    * write output in uncompressed `BAM` format (`-bu`)
* then pipes the output to `msamtools profile` to
    * ignore all multihit inserts (`--multi=ignore`)
    * avoid sequence-length normalization (`--nolen`)
    * calculate actual abundance profiles (`--unit=ab`)
    * use sample label `SAMPLE` in the output file (`--label=SAMPLE`)
    * and write compressed output to `SAMPLE.profile.txt.gz` (`-o`)
