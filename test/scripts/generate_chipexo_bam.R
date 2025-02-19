#!/usr/bin/env Rscript

library(GenomicRanges)
library(Rsamtools)
library(GenomicAlignments)
library(rtracklayer)
library(dplyr)
library(Biostrings)

set.seed(42)  # Ensure reproducibility

# Function to extract Seqinfo from the FASTA file
seqinfo_from_fasta <- function(fasta_file) {
  fasta <- readDNAStringSet(fasta_file)
  seqinfo <- Seqinfo(seqnames = names(fasta), seqlengths = width(fasta))
  return(seqinfo)
}

# Function to extract correct sequence from genome
extract_read_sequence <- function(genome, chr, start, read_length, strand) {
  seq <- subseq(genome[[chr]], start = start, width = read_length)
  
  # Reverse complement for negative strand
  if (strand == "-") {
    seq <- seq
  }
  
  return(as.character(seq))
}

# Function to simulate paired-end reads and generate BAM using GAlignmentPairs
generate_bam <- function(genome, seqinfo_data, annotations, peaks = NULL, output_file, kmer = 5, read_length = 75) {
  fragment_min <- 300
  fragment_max <- 500

  # Total reads per chromosome (uniform coverage)
  total_reads <- 5000  # Read pairs

  read_pairs_first <- list()
  read_pairs_second <- list()

  for (chr in seqlevels(seqinfo_data)) {
    chr_length <- seqlengths(seqinfo_data)[chr]  # Retrieve chromosome length

    # Define a fraction of the chromosome to ignore at both ends
    buffer_size <- floor(0.05 * chr_length)  # 5% of the chromosome

    # Ensure reads are not biased to the start
    r1_start <- sample(seq(buffer_size, chr_length - fragment_max - buffer_size), total_reads, replace = TRUE)

    # Determine fragment sizes
    fragment_sizes <- sample(fragment_min:fragment_max, total_reads, replace = TRUE)

    # Compute second read (R2) start positions
    r2_start <- pmax(r1_start + fragment_sizes - read_length, 1)

    # Compute end positions
    r1_end <- r1_start + read_length - 1
    r2_end <- r2_start + read_length - 1

    # Assign strands
    strand_r1 <- sample(c("+", "-"), total_reads, replace = TRUE)
    strand_r2 <- ifelse(strand_r1 == "+", "-", "+")  # Mate on opposite strand

    # Extract read sequences from genome
    r1_seq <- mapply(extract_read_sequence, MoreArgs = list(genome, chr, read_length),
                     start = r1_start, strand = strand_r1)
    r2_seq <- mapply(extract_read_sequence, MoreArgs = list(genome, chr, read_length),
                     start = r2_start, strand = strand_r2)

    # Create GAlignments objects for R1 and R2
    first_reads <- GAlignments(
      seqnames = Rle(rep(chr, total_reads)),
      pos = r1_start,
      cigar = rep(paste0(read_length, "M"), total_reads),
      strand = strand_r1
    )
    mcols(first_reads)$seq <- r1_seq

    second_reads <- GAlignments(
      seqnames = Rle(rep(chr, total_reads)),
      pos = r2_start,
      cigar = rep(paste0(read_length, "M"), total_reads),
      strand = strand_r2
    )
    mcols(second_reads)$seq <- r2_seq

    read_pairs_first[[chr]] <- first_reads
    read_pairs_second[[chr]] <- second_reads
  }

  # Remove empty elements
  read_pairs_first <- Filter(function(x) length(x) > 0, read_pairs_first)
  read_pairs_second <- Filter(function(x) length(x) > 0, read_pairs_second)

  # Ensure non-empty alignment objects before creating GAlignmentPairs
  if (length(read_pairs_first) == 0 || length(read_pairs_second) == 0) {
    stop("No valid alignments generated!")
  }

  # Merge GAlignments objects safely
  first_all <- Reduce(c, read_pairs_first)
  second_all <- Reduce(c, read_pairs_second)

  # Ensure correct types
  if (!inherits(first_all, "GAlignments") || !inherits(second_all, "GAlignments")) {
    stop("Error: first_all and second_all must be GAlignments objects")
  }

  alignments <- GAlignmentPairs(first_all, second_all)

  # Assign correct seqinfo
  seqinfo(alignments) <- seqinfo_data

  # Create a temporary BAM file for unsorted output
  temp_bam_file <- tempfile(fileext = ".bam")

  # Ensure temporary file is deleted when script exits (even on error)
  on.exit({
    if (file.exists(temp_bam_file)) {
      file.remove(temp_bam_file)
    }
  }, add = TRUE)

  # Export to BAM format
  export(alignments, temp_bam_file, format = "bam")

  # Define final sorted BAM output
  sorted_index_output_name <- paste0(output_file, "_sorted.bam")

  # Sort and index BAM
  sortBam(temp_bam_file, destination=paste0(output_file, "_sorted"))
  indexBam(sorted_index_output_name)

  cat("Generated BAM:", sorted_index_output_name, "\n")
  return(sorted_index_output_name)
}

# Define input genome & annotations
genome_fasta <- "synthetic_genome.fasta"
genome <- readDNAStringSet(genome_fasta)  # Load genome sequences
genome_seqinfo <- seqinfo_from_fasta(genome_fasta)  # Extract reference sequence info

gtf_file <- "synthetic_genome.gtf"

# Load genome & annotations
annotations <- rtracklayer::import(gtf_file)

# Generate Treatment BAM (with peaks at specified genes)
treatment_genes <- c("gene_chr2_1", "gene_chr3_2")  # Example enriched genes
treatment_bam <- generate_bam(genome, genome_seqinfo, annotations, peaks = treatment_genes, output_file = "treatment")

# Generate Control BAM (uniform coverage, no enrichment)
control_bam <- generate_bam(genome, genome_seqinfo, annotations, peaks = NULL, output_file = "control")

cat("BAM files generated:\n")
cat("- Treatment BAM:", treatment_bam, "\n")
cat("- Control BAM:", control_bam, "\n")
