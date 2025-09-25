# Tag Counter

- [Overview](#overview)
- [Input data and Usage](#input-data-and-usage)
- [Installation](#installation)
- [Output](#output)
- [Comparison with Other Tools](#comparison-with-other-tools)
- [Future Enhancements](#future-enhancements)

## Overview

tag_counter is a command-line tool for quantifying the 5′ ends of aligned sequencing
reads over user-defined genomic regions. Unlike traditional peak-calling methods, it
focuses on computing enrichment within explicitly specified regions, such as promoters,
enhancers, or intergenic windows.

## Input data and Usage

1. Genomic Regions (BED6 format)

Input regions should be described in a BED6 file, where each line contains:

<pre>
chr1	463	963	gene1	100	+
chr1	1264	1764	gene2	100	-
chr2	1336	1836	gene3	100	+
chr2	2296	2796	gene4	100	-
chr2	696	1196	gene5	100	+
</pre>

2. 5′ End Coverage Files

These files are produced using bedtools genomecov to count 5′ end tags:  

```bash
bedtools genomecov -ibam <input.bam> -5 -dz > sample_5p_cov.txt
```

Which creates a file like this:

```tsv
chr1	249	1
chr1	250	2
chr1	251	4
chr1	252	1
chr1	256	2
```

Where 

- Column 1: chromosome

- Column 2: 0-based genomic position

- Column 3: number of 5′ ends at that position

These 5' end coverage files are collected into a JSON file for input to `tag_counter`.

```json
{
  "sample1": ["sample1_rep1_5p_cov.txt"],
  "sample2": ["sample2_rep1_5p_cov.txt", "sample2_rep2_5p_cov.txt"]
}
```

If multiple samples can be provided via the json input file. They will be quantified
both individually and and in aggregate.

3. (Optional) Background Sample

If a background sample is provided, then enrichment statistics will be computed. The
background sample would need to be run first in the same way as above, 

```json
{
  "control": ["control_rep1_5p_cov.txt", "control_rep2_5p_cov.txt"]
}
```

and run it first to get the total number of insertions in the background sample.

```bash
tag_counter \
  -b promoters.bed \
  -j control.json \
  -o control
```

But then the total insertions can be found in the file
`<output>/<json_key>/total_tag_counts.tsv` (see Output[#output]). If you provide the
background region counts file (ie `<output>/<json_key>/<json_key>_region_counts.tsv`).
If `tfs.json` is the same as the json described in #2 above, then if you were to run
the following command, then the output for the files in `tfs.json` would contain
enrichment statistics relative to the background sample.

```bash
./tag_counter \
    -b promoters.bed \
    -j tfs.json \
    -o analysis/promoter_quantification \
    --background-counts control/control_region_counts.tsv \
    --background-tag-total 103922
```

## Installation

**tag_counter** is implemented in Rust and can be compiled using Cargo. If you don't
have Rust installed, you can get it from [rustup.rs](https://rustup.rs/).

```bash
git clone https://github.com/YourUsername/tag_counter.git
cd tag_counter
cargo build --release
```

The executable will be located in `target/release/tag_counter`.

## Output

In the specified output directory, there will be subdirectories for each key in the
input JSON file. Each subdirectory will contain the following files:

### `<sample_name>_region_counts.tsv`

A TSV file for each replicate/sample with the following columns:

**Without background:**
- `chrom`: Chromosome name
- `start`: Start position of the region (0-based)
- `end`: End position of the region (1-based)
- `counts`: Number of 5′ end tags in the region

**With background:**
- `chrom`: Chromosome name
- `start`: Start position of the region (0-based)
- `end`: End position of the region (1-based)
- `counts`: Number of 5′ end tags in the region
- `background_counts`: Number of 5′ end tags in the region from the background sample
- `enrichment`: Ratio of normalized counts to normalized background counts
- `p_value`: Poisson p-value testing enrichment significance

### `combined_region_counts.tsv`

Generated only when multiple replicates are provided. Contains the same column structure as individual files but with combined counts across all replicates for each region.

### `total_tag_counts.tsv`

A TSV file with the following columns:

- Column 1: Output filename (e.g., `sample1_region_counts.tsv`) or `combined` for the sum across replicates
- Column 2: Total number of 5′ end tags across the entire genome for that sample
- Column 3: Total number of 5′ end tags within the specified BED regions (unique counts - each tag counted only once even if it overlaps multiple regions)

If multiple replicates are provided, an additional `combined` row shows the totals across all replicates.

## Comparison with Other Tools

### Alternative Approaches for Read Quantification

- One could truncate reads and count using software like
[HTSeq](https://htseq.readthedocs.io/en/latest/counting.html#).
- Truncation may not be necessary **if read lengths are shorter than the minimum
separation between features**.
- If **regions overlap on opposite strands**, truncation is advisable.

### Integration with DESeq2 & EdgeR

- **After quantification**, statistical analysis could be performed using
[DESeq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html) or
[EdgeR](https://bioconductor.org/packages/release/bioc/html/edgeR.html).
- **NOTE**: These packages can be used with the intermediate output of `tag_counter`.

## Planned Enhancements

- replicate comparison metrics
- bootstrapping for confidence intervals of the counts over regions