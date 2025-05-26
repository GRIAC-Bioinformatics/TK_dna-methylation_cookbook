#!/usr/bin/env Rscript
# Author: Vartika Bisht + Tatiana Karp
# Date: May 24, 2025

# Description:
# This script processes a GenomicRatioSet object to remove probes that may be affected by Single Nucleotide Polymorphisms (SNPs).
# SNP information, including Minor Allele Frequency (MAF), is sourced from the dbSNP database, as described in the Bioconductor 450k methylation array documentation.
# The script specifically targets and removes probes that contain a SNP either at the CpG interrogation site or at the single nucleotide extension site.
# An arbitrary decision has been made to filter out probes associated with SNPs having a Minor Allele Frequency (MAF) greater than 10%.

# Usage:
# Rscript 04_snp_containing_flagged_probes.R -g grSet_sample_filtered.Rdata -o low_beadcount_flagged_probes.csv -s 0.1

# Command Line Options:
#   -g, --grset              Path to the genomic ratio set file (.RData).
#   -o, --output             Output file path for the processed GenomicRatioSet (.RData).
#   -s, --snp_maf_threshold  Minor Allele Frequency (MAF) threshold for SNPs (default: 0.1).

suppressPackageStartupMessages({
  tryCatch({
    library(minfi)
    library(dplyr)
    library(optparse)
    library(IlluminaHumanMethylation450kmanifest)
    library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
  }, error = function(e) {
    stop("Package loading failed: ", e$message, call. = FALSE)
  })
})

# Define command line options
option_list <- list(
  make_option(c("-g","--grset"),
              type = "character",
              help = "Path to the genomic ratio set file (.RData)",
              metavar = "PATH"),

  make_option(c("-o","--output"),
              type = "character",
              help = "flagged probes file saved as csv",
              metavar = "FILE"),

  make_option(c("-s","--snp_maf_threshold"),
              type = "numeric",
              default = 0.1,
              help = "Minor Allele Frequency (MAF) threshold for SNPs [default: %default]")
)

# Parse command line arguments
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Use tryCatch for robust error handling
tryCatch({
  # Validate required arguments
  if (is.null(opt$grset) || is.null(opt$output) ) {
    print_help(opt_parser)
    stop("Missing required arguments. Please provide --grset, --output, and --flagged_samples.", call. = FALSE)
  }

  # Load the GenomicRatioSet
  message("Loading GenomicRatioSet from:", opt$grset, "\n")
  load(opt$grset) # This should load an object named 'grSet' into the environment

  # Check if grSet object is loaded
  if (!exists("grSet") || !inherits(grSet, "GenomicRatioSet")) {
    stop("The loaded .RData file does not contain a 'GenomicRatioSet' object named 'grSet'.", call. = FALSE)
  }

  # Store all initial probe names before dropping
  all_initial_probes <- rownames(grSet)
  message("Initial number of probes:", length(all_initial_probes), "\n")

  # Drop probes containing SNPs
  message("Dropping probes with SNPs based on MAF threshold >", opt$snp_maf_threshold, "\n")
  grSet_filtered <- dropLociWithSnps(grSet, snps = c("SBE", "CpG"), maf = opt$snp_maf_threshold)

  # Identify flagged probes
  removed_probes <- all_initial_probes[!(all_initial_probes %in% rownames(grSet_filtered))]
  message("Number of probes removed due to SNPs:", length(removed_probes), "\n")
  message("Number of probes remaining after SNP filtering:", nrow(grSet_filtered), "\n")

  # Save the flagged probes
  message("Saving flagged probes to :", opt$output, "\n")
  write.csv(data.frame(rownames(grSet_filtered)), file = opt$output, row.names = FALSE, 
            col.names = "Flagged_Probes") 

}, error = function(e) {
  # messagech and report any errors
  message("An error occurred during script execution:\n")
  message(e$message, "\n")
  stop("Script terminated due to error.", call. = FALSE)
})