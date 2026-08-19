# msamtools:  microbiome-related extension to samtools

**msamtools** provides useful functions that are commonly used in microbiome
data analysis, especially when analyzing shotgun metagenomics
or metatranscriptomics data.

# Using msamtools

This is the master program that you call with the subprogram options. There
are currently 4 subprograms that you can call as shown below.

```text
Program: msamtools (Metagenomics-related extension to samtools)
Version: 1.1.3 (git 3fd8379a7189; using htslib 1.24)

Usage:   msamtools <command> [options]

Commands:
 -- Filtering
     filter         filter alignments based on alignment statistics

 -- Profiling
     profile        estimate relative abundance profile of reference sequences or genomes in bam file

 -- Coverage
     coverage       estimate per-base or per-sequence read coverage of each reference sequence

 -- Summary
     summary        summarize alignment statistics per read in a table format
```

These represent the different analysis of SAM/BAM files you can perform using
**msamtools**.
The [example workflows](examples.md) show how to combine them with each other or **samtools**. The command pages explain [filter](filter.md), [profile](profile.md), [coverage](coverage.md), and [summary](summary.md) in detail.
