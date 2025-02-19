#!/usr/bin/env Rscript

library(Biostrings)
library(GenomicRanges)
library(GenomicFeatures)
library(rtracklayer)
library(Rsamtools)
library(dplyr)

set.seed(42)  # Ensure reproducibility

# Function to generate a random chromosome sequence
generate_chromosome <- function(length) {
  bases <- c("A", "C", "G", "T")
  DNAString(paste(sample(bases, length, replace = TRUE), collapse = ""))
}

# Function to generate gene annotations with placement restrictions
generate_gene_annotations <- function(chr_name, chr_length, gene_length = 300) {
  min_pos <- 800  # Ensure genes start at least 800 bp from the start
  max_pos <- chr_length - 800 - gene_length  # Ensure genes end at least 800 bp from the end

  # First gene: placed randomly within valid range
  pos1 <- sample(min_pos:max_pos, 1)

  # Second gene: overlapping on opposite strand
  overlap_start <- sample((pos1 - 100):(pos1 + gene_length - 100), 1)
  overlap_start <- max(min_pos, min(overlap_start, max_pos))

  genes <- GRanges(
    seqnames = Rle(chr_name, 2),
    ranges = IRanges(start = c(pos1, overlap_start),
                     end = c(pos1 + gene_length - 1, overlap_start + gene_length - 1)),
    strand = c("+", "-"),
    source = "synthetic",
    feature = "gene",
    score = NA,
    frame = NA,
    gene_id = paste0("gene_", chr_name, "_", c(1, 2))
  )

  # Chromosome 2: Add a third feature upstream of gene_chr2_2, on opposite strand
  if (chr_name == "chr2") {
    new_gene_start <- max(min_pos, start(genes[2]) - 500 - gene_length)  # 500bp upstream of gene_chr2_2
    new_gene_end <- new_gene_start + gene_length - 1
    new_gene <- GRanges(
      seqnames = chr_name,
      ranges = IRanges(start = new_gene_start, end = new_gene_end),
      strand = ifelse(strand(genes[2]) == "+", "-", "+"),  # Opposite strand
      source = "synthetic",
      feature = "gene",
      score = NA,
      frame = NA,
      gene_id = paste0("gene_", chr_name, "_3")
    )
    genes <- c(genes, new_gene)
  }

  # Chromosome 3: Same logic as chromosome 2, but upstream of gene_chr3_1
  if (chr_name == "chr3") {
    new_gene_start <- max(min_pos, start(genes[1]) - 500 - gene_length)  # 500bp upstream of gene_chr3_1
    new_gene_end <- new_gene_start + gene_length - 1
    new_gene <- GRanges(
      seqnames = chr_name,
      ranges = IRanges(start = new_gene_start, end = new_gene_end),
      strand = ifelse(strand(genes[1]) == "+", "-", "+"),  # Opposite strand
      source = "synthetic",
      feature = "gene",
      score = NA,
      frame = NA,
      gene_id = paste0("gene_", chr_name, "_3")
    )
    genes <- c(genes, new_gene)
  }

  list(genes = genes)
}

# Define genome structure (chromosomes of 5000 bp)
chromosome_lengths <- c(chr1 = 5000, chr2 = 5000, chr3 = 5000)

# Generate synthetic genome sequences
synthetic_genome <- lapply(names(chromosome_lengths), function(chr) {
  generate_chromosome(chromosome_lengths[chr])
})
names(synthetic_genome) <- names(chromosome_lengths)

# Generate GRanges object for genes
annotations <- lapply(names(synthetic_genome), function(chr) {
  generate_gene_annotations(chr, chromosome_lengths[chr])
})

gene_annotations <- do.call(c, lapply(annotations, `[[`, "genes"))

# Ensure seqlevels match between genome and annotations
seqlevels(gene_annotations) <- names(synthetic_genome)

# Write genome sequences to a FASTA file
genome_fasta_file <- "synthetic_genome.fasta"
writeXStringSet(DNAStringSet(synthetic_genome), filepath = genome_fasta_file, format = "fasta")

# Create FASTA index (.fai) file
indexFa(genome_fasta_file)

# File paths for annotations
gtf_file <- "synthetic_genome.gtf"
gff3_file <- "synthetic_genome.gff3"
bed_file <- "synthetic_genome_promoters.bed"

# Generate promoter regions using GenomicFeatures::promoters()
promoters <- promoters(gene_annotations, upstream = 500, downstream = 0)

# Ensure promoters are within valid regions (not before position 1)
promoters <- promoters[start(promoters) >= 1]

# Assign metadata for BED format
mcols(promoters)$score <- 100  # Arbitrary score
mcols(promoters)$gene_id <- mcols(gene_annotations)$gene_id  # Assign gene ID

# Ensure seqlevels are consistent
seqlevels(promoters) <- seqlevels(gene_annotations)

# Export GTF and GFF3 files
export(gene_annotations, con = gtf_file, format = "gtf")
export(gene_annotations, con = gff3_file, format = "gff3")

# Convert promoters to a BED-compatible data frame
promoter_bed <- data.frame(
  chrom = seqnames(promoters),
  start = start(promoters) - 1,  # BED is 0-based
  end = end(promoters),
  name = mcols(promoters)$gene_id,
  score = mcols(promoters)$score,
  strand = strand(promoters)
)

# Write BED file
write.table(promoter_bed, bed_file, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

# Print output files
cat("Synthetic genome and annotation files created:\n")
cat("- Genome file:", genome_fasta_file, "\n")
cat("- FASTA index file:", paste0(genome_fasta_file, ".fai"), "\n")
cat("- GTF annotation file:", gtf_file, "\n")
cat("- GFF3 annotation file:", gff3_file, "\n")
cat("- BED promoter file:", bed_file, "\n")

