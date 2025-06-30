#!/usr/bin/env Rscript
# Script: 01_create_sigset_object.R
# Author: Vartika Bisht
# Date: June 22, 2025
#
# Purpose: Process Illumina methylation array IDAT files into RGChannelSet
#
# Input: Directory containing IDAT files
# Output: RGChannelSet object saved as RData file
#
# Usage:
#    Rscript 01_create_sigset_object.R --basedir /path/to/idats --output output.RData --platform HM450
#


# This script processes IDAT files using the 'sesame' R package.
# It expects two command-line arguments:
# 1. Base directory containing IDAT files.
# 2. Output file path for the generated SigSet object (RData format).

# | File HM450.address needs to be cached to be used in sesame.
# | Please make sure you have updated ExperimentHub and try
# | > sesameDataCache()
# | to retrieve and cache needed sesame data.
# 
# https://rdrr.io/github/zwdzwd/sesame/f/vignettes/sesame.Rmd

# install.packages("devtools") # If you don't have devtools
# devtools::install_version("dbplyr", version = "2.3.4")

# Load necessary libraries

suppressPackageStartupMessages({
  tryCatch({
  library(sesame)
  library(optparse)
  library(BiocParallel)
}, error = function(e) {
    stop("Package loading failed: ", e$message, call. = FALSE)
  })
})

# You need to cache the sesame datasets first before you start.
# This has to be done only once. If you have already done this, you can skip this step.
# In bash declare the following variables:
# export EXPERIMENT_HUB_CACHE="/groups/umcg-griac/tmp02/projects/Vartika/R_Cache/"
# export BIOCFILECACHE="/groups/umcg-griac/tmp02/projects/Vartika/R_Cache/"
# Then in R run the following commands:
# library(ExperimentHub)
# library(BiocFileCache)
# sesameDataCache()


# Define command-line options
option_list = list(
  make_option(c("-d", "--basedir"), type="character", default=NULL, 
              help="base directory containing IDAT files (e.g., /groups/umcg-griac/tmp02/rawdata/medall/IDAT)", metavar="character"),
  make_option(c("-p", "--platform"), type="character", default=NULL, 
              help="platform type (e.g., HM450)", metavar="character"),
  make_option(c("-c", "--cores"), type="character", default=4, 
              help="number of cores to use (e.g., 4)", metavar="character"),
  make_option(c("-o", "--output_sigset_file"), type="character", default=NULL, 
              help="output file name for SigSet object (e.g., output.Rdata)", metavar="character")
); 

# Parse command-line arguments
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

# Check if required arguments are provided
if (is.null(opt$basedir)){
  print_help(opt_parser)
  stop("Base directory argument must be supplied (-d/--basedir)", call.=FALSE)
}
if (is.null(opt$output_sigset_file)){
  print_help(opt_parser)
  stop("Output SigSet file argument must be supplied (-o/--output_sigset_file)", call.=FALSE)
}

# Assign base directory from parsed options
dir.name <- opt$basedir

# Use tryCatch for robust error handling during file operations
tryCatch({
  # Search for IDAT file prefixes within the specified directory
  message(paste("Searching for IDAT prefixes in:", dir.name))
  idats <- searchIDATprefixes(dir.name, recursive = TRUE, use.basename = TRUE)
  
  # Check if any IDAT files were found
  if (length(idats) == 0) {
    stop("No IDAT files found in the specified base directory: ", dir.name)
  }
  
  # Determine the number of cores to use for parallel processing
  message(paste("Using", opt$cores, "cores for parallel IDAT reading."))
  
  # Set up the parallel backend. MulticoreParam is suitable for Unix-like systems.
  # For Windows, SnowParam or registerDoSNOW.
  register(MulticoreParam(opt$cores))
  
  # Read each IDAT pair into a SigSet object using BiocParallel for parallel processing.
  # This creates a list where each element is a SigSet object corresponding to one IDAT pair.
  message(paste("Reading all IDAT pairs into SigSet objects for platform:", opt$platform))
  SigSetList <- bplapply(idats, function(pfx) {
    message("Processing sample:", pfx$sampleName)
    readIDATpair(pfx, platform = opt$platform, manifest = NULL, controls = NULL, verbose = FALSE)
  }, BPPARAM = bpparam()) # Use the registered parallel backend

  # Save the SigSet object to the specified output file
  message(paste("Saving SigSet object to:", opt$output_sigset_file))
  save(SigSetList, file = opt$output_sigset_file)
  message("SigSet object saved successfully.")
}, error = function(e) {    
  # Handle errors and stop execution with an informative message
  stop("Failed to process IDAT files or save output file: ", e$message)
})
