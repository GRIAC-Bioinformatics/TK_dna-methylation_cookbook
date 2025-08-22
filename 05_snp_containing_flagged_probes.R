#!/usr/bin/env Rscript
# =============================================================================
# 05_snp_containing_flagged_probes.R
# Author: Vartika Bisht & Tatiana Karp
# Created: 13-Aug-2025
# Description:
# This script filters out probes that may have compromised signal integrity due
# to the presence of a Single Nucleotide Polymorphism (SNP). It specifically
# removes probes where a known SNP (with a Minor Allele Frequency > 10%) is
# located at the CpG interrogation site or the single nucleotide extension site.
# Arguments:
#    --grset                → Path to the GenomicRatioSet object file (.RData)
#    --flagged              → Output CSV file to list flagged probes
#    --cutoff               → Minor Allele Frequency (MAF) threshold for SNPs (default: 0.1)
# Usage:
#   Rscript 05_snp_containing_flagged_probes.R --grset <grSet.RData> \
#     --flagged <snp_flagged_probes.csv> --cutoff 0.1
# Notes:
#   - The presence of a SNP can disrupt probe hybridization and/or the single-base
#     extension, leading to inaccurate methylation measurements.
#   - This script uses SNP information, typically from databases like dbSNP, to
#     identify problematic probes.
#   - The filtering is based on the Minor Allele Frequency (MAF). A higher
#     MAF indicates that the SNP is more common in the population, making it more
#     likely to affect your samples. The default cutoff of 0.1 (10%) is a
#     common, though arbitrary, choice.
#   - Filtering these probes is a critical step to ensure that methylation
#     calls are not confounded by genetic variation.
# =============================================================================

suppressPackageStartupMessages({
  tryCatch({
    library(minfi)
    library(dplyr)
    library(optparse)
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

  make_option(c("-f","--flagged"),
              type = "character",
              help = "flagged probes file saved as csv",
              metavar = "FILE"),

  make_option(c("-c","--cutoff"),
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
  if (is.null(opt$grset) || is.null(opt$flagged) ) {
    print_help(opt_parser)
    stop("Missing required arguments. Please provide --grset, --flagged, and --flagged_samples.", call. = FALSE)
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
  message("Dropping probes with SNPs based on MAF threshold >", opt$cutoff, "\n")
  grSet_filtered <- dropLociWithSnps(grSet, snps = c("SBE", "CpG"), maf = opt$cutoff)

  # Identify flagged probes
  removed_probes <- all_initial_probes[!(all_initial_probes %in% rownames(grSet_filtered))]
  message("Number of probes removed due to SNPs:", length(removed_probes), "\n")
  message("Number of probes remaining after SNP filtering:", nrow(grSet_filtered), "\n")

  # Save the flagged probes
  message("Saving flagged probes to :", opt$flagged, "\n")
  write.csv(data.frame("Flagged_Probes" = removed_probes), file = opt$flagged, row.names = FALSE) 

}, error = function(e) {
  # messagech and report any errors
  message("An error occurred during script execution:\n")
  message(e$message, "\n")
  stop("Script terminated due to error.", call. = FALSE)
})