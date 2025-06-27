# Tag Counter

- [Overview](#overview)
- [Installation](#installation)
- [Usage](#usage)
- [Algorithmic Details](#algorithmic-details)
- [Comparison with Other Tools](#comparison-with-other-tools)
- [Future Enhancements](#future-enhancements)

## Overview

**tag_counter** is a tool for analyzing genomic sequencing data that generates
positional tags at the 5' end of aligned reads. Unlike peak-calling
methods, **tag_counter** focuses on enrichment over user defined genomic regions.
This tool quantifies enrichment, e.g. of reads produced in a
[callingcards](https://nf-co.re/blog/2024/callingcards_intro) experiment for a given
transcription factor, by comparing treatment and control samples to produce statistical
enrichment values and quality control (QC) metrics.

## Installation

**tag_counter** is implemented in Rust and can be compiled using Cargo.

```bash
cargo build --release
```

## Usage

To run **tag_counter**:

```bash
tag_counter \
  --bed-file synthetic_genome_promoters.bed \
  --json-file input.json
  --regions promoters.bed \
  --output example_tag_counter_output 
```

### Required Arguments

- `--treatment`: One or more BAM/CRAM files containing the treatment reads.
- `--control`: One or more BAM/CRAM files containing the control reads.
- `--regions`: A BED file containing the genomic regions to analyze.

### Optional Arguments
- `--bootstrap-sample-size`: The number of reads to sample for each bootstrap sample.
(default: 1e6)
- `--unstratified-bootstrap`: Setting this flag disables the stratified bootstrap
sampling. See [Generate Bootstrap Samples](#generate-bootstrap-samples) for more
details. Stratified sampling is enabled by default.
- `--kmer`: The length of the k-mer (tag length from the 5' end) to analyze
(default: 5).
- `--output`: The name of the output file to save the results.

## Algorithmic Details

**tag_counter** assumes that the input BAM/CRAM files contain preprocessed and
filtered reads, meaning any prior filtering steps (e.g., quality filtering,
duplicate removal) should already be completed before running the tool.

### Generate Bootstrap Samples

First, calculate the total number of reads in each replicate. Then, determine
how many reads to sample from each replicate based on the bootstrap sample size. For 
example, if there is only one file, then 100% of the reads will originate from that
file. If there are two files, and the first replicate has 15e6 reads while the
second has 5e6 reads, then the first replicate will contribute 75% of the reads to the
bootstrap sample, while the second will contribute 25%. These values are rounded so
that the total number of reads in each bootstrap samples is equal to the value set by
`--bootstrap-sample-size`.

Within each of the replicate files, reads are sampled in a stratified manner, with
replacement in the following manner:

- Reads are first categorized into region-associated and non-region-associated groups.
- Within each category, reads are further stratified based on coverage quantiles
(low-, medium-, and high-coverage regions).
- Bootstrap resampling is performed proportionally from each coverage tier,
ensuring that low-coverage regions contribute meaningfully to enrichment
variability estimates.
- This approach prevents highly covered regions from dominating the bootstrap
resampling, leading to more reliable and unbiased confidence intervals for
enrichment estimates.

### Region and K-mer Quantification

- Iterate over treatment reads, and do the following:
    - If the first k nucleotides overlap one of the user specified
    regions, then that region's count is incremented.
    - Record the composition of the first k nucleotides for k-mer analysis as follows:
        - Increment the k-mer count for the 'global' k-mer structure
        - if the k-mer overlaps a region, increment the count of the 'region' k-mer
        structure
        - if the k-mer does not overlap a region, increment the count of the
        'non-region' k-mer structure
- In parallel, if there are adequate resources, the same process is performed on
the control reads.
- In parallel, if there are adequate resources, the same process is performed on
the treatment and control bootstrap samples.

### Enrichment Analysis

For each region, tag_counter calculates the enrichment score, Poisson p-value,
hypergeometric p-value, and a confidence interval derived from bootstrap
samples.

### K-mer Analysis

The nucleotide frequency of the k-mers is calculated for the treatment and control for
each of the three categories: global, region, and non-region. The variance of the 
k-mer is estimated from the bootstrap samples. A high variance suggests potential
biases in replicate composition. If no replicates are present, high variance may
indicate a bias in read composition between treatment and control. The categories are
also compared to one another, with the expectation that in the control, the three would
be similar, while in the treatment we would expect that region vs non-region to be
different. Finally, the categories are compared between treatment and control. We
expect that the region k-mers in the treatment will be different from the region
k-mers in the control. If the non-region k-mers are different between treatment and
control, and the region and non-region k-mers are different from one another within the
treatment, then further investigation is warranted to determine if there is a bias in
the treatment sample that isn't accounted for in the control.

## Comparison with Other Tools

A similar analysis could be performed by **k-mer analysis of the 5' ends** of
the filtered reads. An unfiltered analysis would provide **more insight into
sequencing and/or library prep artifacts**. However, if the goal is to detect
**artifacts that remain after filtering**, then analyzing k-mer content on the
**filtered reads** is appropriate.

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

- **Python and R bindings** for easier integration into bioinformatics workflows.
- **Peak calling methods** for different binding assay modalities. This would
provide the ability to do more traditional peak calling and compare the overlap
between peaks and the user provided regions
- **A Control free method** which would estimate the noise/bias from the treatment
files only
- **Internal implementation of an interval tree** to remove the `bio` dependency