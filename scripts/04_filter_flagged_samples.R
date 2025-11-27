#!/usr/bin/env Rscript
# =============================================================================
# 04_filter_flagged_samples.R
# Author: Vartika Bisht
# Created: 13-Aug-2025
# Description:
# This script filters samples from one or more minfi R data objects based on
# a list of flagged samples from a quality control analysis. It removes samples
# that have been flagged two or more times, saving the cleaned objects to new
# files with a user-specified suffix.
# Arguments:
#    --flagged              → Path to a CSV file containing flagged samples (e.g., from 04_sample_qc_overview.R)
#    --mset                 → Path to the MSet.RData file (optional)
#    --grset                → Path to the GenomicRatioSet.RData file (optional)
#    --rgset                → Path to the RGChannelSet.RData file (optional)
#    --rgsetext             → Path to the RGChannelSetExtended.RData file (optional)
#    --base_suffix          → Suffix to add to the output file names (default: _sample_filtered)
#    --pdf_output           → Output PDF file path for plots (optional)
# Usage:
#   Rscript 04_filter_flagged_samples.R --flagged <flagged_samples.csv> \
#     --mset <MSet.RData> --rgset <RGset.RData> --base_suffix "_filtered_samples"
# Notes:
#   - The script reads a CSV file containing a list of samples and their total number of flags.
#     Samples with a Total_Occurrences of 2 or more are identified as "bad" and are removed.
#   - This is a crucial step for data cleaning after initial QC, ensuring that subsequent
#     analyses are performed on high-quality samples.
#   - The output files are saved in the same directory as their input counterparts.
# =============================================================================

suppressPackageStartupMessages({
  tryCatch({
library(optparse)
library(minfi) 
}, error = function(e) {
    stop("Package loading failed: ", e$message, call. = FALSE)
  })
})

# Parse command line arguments
option_list <- list(
  make_option(c("-s", "--flagged"), type = "character", default = NULL,
              help = "Path to the CSV file with list of samples to be excluded."),
  make_option(c("-m", "--mset"), type = "character", default = NULL,
              help = "Path to the MSet.RData file"),
  make_option(c("-g", "--grset"), type = "character", default = NULL,
              help = "Path to the grSet.RData file"),
  make_option(c("-R", "--rgset"), type = "character", default = NULL,
              help = "Path to the RGset.RData file"),
  make_option(c("-x", "--rgsetext"), type = "character", default = NULL,
              help = "Path to the RGsetEXT.RData file"),
  make_option(c("-b", "--base_suffix"), type = "character", default = "_sample_filtered",
              help = "Suffix to add to output file names (e.g., '_filtered'). Default is '_sample_filtered'."),
  make_option(c("-p", "--pdf_output"), type = "character", default = NULL,
              help = "Name for the output PDF plot file containing all plots (required).")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# VALIDATE REQUIRED ARGUMENTS
required_args <- c("mset", "grset", "rgset", "rgsetext" , "base_suffix" , "flagged" , "pdf_output")

for (arg in required_args) {
  if (is.null(opt[[arg]])) {
    stop(paste("Error: Required argument --", arg, " is missing. Use --help for more information.", sep = ""))
  }
}


# Collect all R data file paths into a named list for easier iteration
# The name in the list should match the object name loaded from the .RData file
input_rdata_paths <- list(
  MSet = opt$mset,
  grSet = opt$grset,
  RGset = opt$rgset,
  RGsetEXT = opt$rgsetext
)

# Filter out NULL entries for any R data files not specified by the user
input_rdata_paths <- input_rdata_paths[!sapply(input_rdata_paths, is.null)]

# Check if all specified R data files actually exist on disk
for (file_path in input_rdata_paths) {
  if (!file.exists(file_path)) {
    stop(paste("The specified R data file does not exist:", file_path))
  }
}

### Identify Samples to Exclude

# Read the sample QC CSV file into a data frame
flagged <- read.csv(opt$flagged, header = FALSE)[[1]]

message(length(flagged), " samples to be excluded.")

### Filter Each R Data Object

process_minfi_object <- function(file_path, object_name, samples_to_exclude, base_suffix) {
  message(paste("\nProcessing file:", file_path))

  # Load the R data object. It will create a variable with 'object_name' (e.g., 'MSet').
  load(file_path)

  # Get the loaded object dynamically by its name.
  current_object <- get(object_name)

  # Determine the indices of samples to keep.
  samples_to_keep_idx <- which(!sampleNames(current_object) %in% samples_to_exclude)

  # Subset the object to keep only the desired samples.
  # Assign it directly back to the original object_name variable.
  assign(object_name, current_object[, samples_to_keep_idx])

  # Construct the output file name.
  base_file_name <- tools::file_path_sans_ext(basename(file_path))
  output_file_name <- file.path(dirname(file_path), paste0(base_file_name, base_suffix, ".RData"))

  # Save the *modified* object (which has the original 'object_name' variable name).
  save(list = object_name, file = output_file_name)
  message(paste("  Filtered and saved:", output_file_name))
}

# Process MSet
process_minfi_object(opt$mset, "MSet", flagged, opt$base_suffix)

# Process grSet
process_minfi_object(opt$grset, "grSet", flagged, opt$base_suffix)

# Process RGset
process_minfi_object(opt$rgset, "RGset", flagged, opt$base_suffix)

# Process RGsetEXT
process_minfi_object(opt$rgsetext, "RGsetEXT", flagged, opt$base_suffix)


message("\nSample filtering process completed for all specified R data files.")

## Overall distributions of Beta values for each sample
# Here we expect beta values presenting values close to zero or one. 
# This is an indication of unmethylatilated CpG sites when close to zero and methylated when close to one.
# if you have NAs due to 0 intesity values you will have errors here so to avoid this, use offset

tryCatch({
    message("Plotting overall distributions of Beta values for each sample..")

    # Open PDF device for plotting
    pdf(opt$pdf_output, width = 8, height = 6)
    message(paste("PDF plot will be saved to:", opt$pdf_output))

    # --- Load the initial MSet object ---
    message(paste("Loading initial MSet object from:", opt$mset))
    
    load(opt$mset)

    if (!exists("MSet")) {
        stop("Error: 'MSet' object not found after loading initial MSet file.")
    }

    # Extract beta values for the original MSet
    beta_original_matrix <- getBeta(MSet, offset = 100)
    message("Successfully extracted beta values from initial MSet.")
    
    rm(MSet)
    gc()

    # Plot density plot
    print(minfi::densityPlot(
      dat = beta_original_matrix,
      main = "Density Plot of Beta Values:Raw",
      xlab = "Beta Value",
      xlim = c(0, 1),
      legend = TRUE
    ))

    rm(beta_original_matrix)

    # --- Load the filtered_object (which will be renamed to MSet) ---
    filtered_mset_path <- file.path(dirname(opt$mset), paste0(tools::file_path_sans_ext(basename(opt$mset)), opt$base_suffix, ".RData"))

    message(paste("Loading filtered MSet object from:", filtered_mset_path))

    load(filtered_mset_path) 

    # Extract beta values for the filtered MSet
    beta_filtered_matrix <- getBeta(MSet, offset = 100)
    message("Successfully extracted beta values from filtered MSet.")

    # Remove the filtered MSet object as its beta values are now extracted
    rm(MSet)

    # Plot density plot
    print(minfi::densityPlot(
      dat = beta_filtered_matrix,
      main = "Density Plot of Beta Values:Filtered",
      xlab = "Beta Value",
      xlim = c(0, 1),
      legend = TRUE
    ))

    dev.off()
}, error = function(e) {
    message("An error occurred during plotting: ", e$message)
})

message("Plot saved.")