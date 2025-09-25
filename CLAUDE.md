# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

`tag_counter` is a Rust command-line tool for quantifying 5′ ends of aligned sequencing reads over user-defined genomic regions. It focuses on computing enrichment within explicitly specified regions (promoters, enhancers, etc.) rather than traditional peak-calling methods.

The tool processes:
- BED6 format genomic regions
- 5′ end coverage files (from `bedtools genomecov`)
- Optional background samples for enrichment statistics

## Commands

### Build and Development
```bash
# Build the project
cargo build --release

# Run tests
cargo test --verbose

# Build with debug info
cargo build --verbose
```

### Running the Tool
```bash
# Basic usage
./target/release/tag_counter -b regions.bed -j samples.json -o output_dir

# With background enrichment
./tag_counter \
    -b promoters.bed \
    -j treatment.json \
    -o analysis/promoter_quantification \
    --background-counts control/control_region_counts.tsv \
    --background-tag-total 103922
```

### Test Data Generation
Test data generation scripts are available in `tests/scripts/`:
- `generate_chipexo_bam.R` - Generate ChIP-exo BAM files
- `generate_fasta_gtf.R` - Generate FASTA and GTF files

## Code Architecture

### Core Components

**`src/main.rs`** - Command-line interface and main processing logic
- Argument parsing with clap
- JSON input file processing
- Background enrichment calculations
- Output file generation

**`src/lib.rs`** - Core data structures and algorithms
- `CountableRegionTree`: Main data structure using interval trees for efficient genomic region queries
- `GenomicInterval<T>`: Represents genomic intervals with optional metadata
- `CountableRegionMetadata`: Stores replicate counts and coverage levels

### Module Structure

**`src/parsers/`** - Input file parsing
- `parse_bed6.rs` - BED6 format parser
- `tag_count_parser.rs` - 5′ end coverage file parser
- `region_count_parser.rs` - Background region count file parser

**`src/quantification/`** - Statistical analysis
- `region_stats.rs` - Enrichment and Poisson p-value calculations

**`src/intervaltree/`** - Interval tree implementation for efficient overlap queries
- Custom interval tree with genomic coordinate support
- Non-overlapping interval construction from overlapping input regions

### Key Algorithms

**Interval Tree Construction**: The tool converts potentially overlapping user-defined regions into non-overlapping canonical intervals using BTreeSet for efficient querying.

**Count Aggregation**: Uses interval trees to efficiently map 5′ end positions to genomic regions, handling multiple replicates simultaneously.

**Statistical Analysis**: Computes enrichment ratios and Poisson p-values when background data is provided.

## Dependencies

- **clap**: Command-line argument parsing
- **serde/serde_json**: JSON input file parsing
- **statrs**: Statistical calculations (Poisson distribution)
- **thiserror**: Error handling
- **pyo3**: Python bindings (optional)
- **extendr-api**: R bindings (optional)

## Input/Output Formats

**Input**:
- BED6 files: `chr1 463 963 gene1 100 +`
- JSON sample mappings: `{"sample1": ["file1_5p_cov.txt"]}`
- 5′ coverage files: `chr1 249 1` (chrom, 0-based position, count)

**Output**:
- `<sample>_region_counts.tsv`: Per-region counts with optional enrichment
- `total_tag_counts.tsv`: Total counts per sample
- `combined_region_counts.tsv`: Aggregated counts (when multiple replicates)

## Testing

Tests are located in `src/lib.rs` and cover:
- Interval tree construction and querying
- Count aggregation across replicates
- Edge cases (no overlaps, multiple overlaps)
- Statistical calculations