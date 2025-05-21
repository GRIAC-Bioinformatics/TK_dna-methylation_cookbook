#!/usr/bin/env Rscript
###############################################################################
# Title: Methylation Array Data Preparation Pipeline
# Author: Vartika Bisht
# Date: 20 May 2025
#
# Purpose: Process Illumina methylation array IDAT files into RGChannelSet
#
# Input: Directory containing IDAT files and sample sheet
# Output: RGChannelSet object saved as RData file
#
# Usage:
#   Rscript methylation_preprocessing.R --basedir /path/to/idats --output output.RData
#
# Memory Requirements:
#   - Expect ~20GB RAM for 450k array with 22 samples
#   - Parallel processing not recommended due to memory constraints
###############################################################################

# Load required packages silently with error handling
suppressPackageStartupMessages({
  required_packages <- c(
    "minfi",
    "optparse",
    "IlluminaHumanMethylation450kmanifest",
    "IlluminaHumanMethylation450kanno.ilmn12.hg19"
  )
  
  for (pkg in required_packages) {
    if (!require(pkg, character.only = TRUE)) {
      stop("Package ", pkg, " not available. Please install from Bioconductor.")
    }
  }
})


option_list <- list(
  make_option(c("--basedir"), 
              type = "character",
              help = "Base directory containing IDAT files and sample sheet",
              metavar = "PATH"),
  
  make_option(c("--output"), 
              type = "character",
              help = "Output file path for RGChannelSet (.RData)",
              metavar = "FILE")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Validate required arguments
if (is.null(opt$basedir) || is.null(opt$output)) {
  print_help(opt_parser)
  stop("Both --basedir and --output arguments must be provided.")
}


# SAMPLE SHEET PROCESSING
message("Reading sample sheet...")

targets <- tryCatch({
  minfi::read.metharray.sheet(
    base = opt$basedir,
    pattern = "Samplesheet.*\\.csv$",
    recursive = TRUE,
    verbose = TRUE
  )
}, error = function(e) {
  stop("Failed to read sample sheet: ", e$message)
})

message(sprintf(
  "Successfully loaded sample sheet with %d samples", 
  nrow(targets)
))

# IDAT FILE PROCESSING
message("Reading IDAT files (this may take several minutes)...")

# read.metharray.exp() fully loads all IDATs into memory and constructs dense matrices.
# This is memory-intensive due to duplicated quant matrices and large final RGChannelSet.
# Expect RAM usage to exceed 2x the IDAT file size. For 22 samples (450k), >60GB is normal.
# List of Quant matrices → Combined Matrix:
# GreenMean <- do.call(cbind, lapply(G.Quants, ...))
# RedMean   <- do.call(cbind, lapply(R.Quants, ...))
# This constructs large in-memory double-precision numeric matrices.
# These matrices are dense and not compressed. For ~485,000 probes × 22 samples × 2 channels (Green + Red), this is a lot of RAM.
# https://rdrr.io/bioc/minfi/man/read.metharray.exp.html 
####################################################################################################################################################
# Instead of using read.metharray.exp(), we will use read.metharray.exp.par() which constructs RGChannelSet in parallel using BiocParallel.
# BiocParallel uses MulticoreParam on Unix-like systems (Linux/macOS), which distributes the workload across multiple CPU cores on a single machine.
# The function splits the IDAT files into chunks and processes them in parallel. Each core reads and parses a subset of the IDAT files independently. 
# After all files are read, the results are combined into a single RGChannelSet object.
####################################################################################################################################################
# Well, that did not work. the par version fo read.metharray.exp() always runs out of memory. If the files are too large, it will crash as you need
# quite some RAM to read the IDAT files. So , you do not want to do it parallel.

RGset <- tryCatch({
  read.metharray.exp(targets = targets)
}, error = function(e) {
  stop("IDAT processing failed. Possible memory issues? Error: ", e$message)
})

# SAVE OUTPUT
message("\nSaving RGChannelSet object...")

tryCatch({
  save(RGset, file = opt$output)
}, error = function(e) {
  stop("Failed to save output file: ", e$message)
})
