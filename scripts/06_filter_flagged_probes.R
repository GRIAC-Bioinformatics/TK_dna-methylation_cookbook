#!/usr/bin/env Rscript
# =============================================================================
# 06_filter_flagged_probes.R
# Author: Vartika Bisht
# Created: 13-Aug-2025
# Description:
# This script filters probes from one or more minfi R data objects based on
# a list of flagged probes from a quality control analysis. It removes probes
# that have been flagged for various reasons, saving the cleaned objects to new
# files with a user-specified suffix. This is a critical step for data cleaning
# before normalization and downstream analysis.
# Arguments:
#    --flagged              → Path to a CSV file containing the list of probes to be excluded
#    --mset                 → Path to the MSet.RData file (optional)
#    --grset                → Path to the GenomicRatioSet.RData file (optional)
#    --rgset                → Path to the RGChannelSet.RData file (optional)
#    --rgsetext             → Path to the RGChannelSetExtended.RData file (optional)
#    --base_suffix          → Suffix to add to the output file names (default: _probe_filtered)
#    --pdf_output           → Output PDF file path for plots (optional)
# Usage:
#   Rscript 06_filter_flagged_probes.R --flagged <flagged_probes.csv> \
#     --mset <MSet.RData> --grset <grSet.RData> --base_suffix "_filtered_probes"
# Notes:
#   - This script takes the output from probe-level QC scripts (e.g., those flagging
#     low bead count, high intensity, or SNP-containing probes) and removes the
#     corresponding probes from the specified data objects.
#   - The _probe_filtered suffix is added to the output file names to distinguish
#     them from the original files.
#   - This step is essential for ensuring that downstream normalization and statistical
#     analyses are performed only on high-quality, reliable probes.
# =============================================================================

suppressPackageStartupMessages({
  tryCatch({
    library(minfi)
    library(optparse)
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

### Identify probes to Exclude

# Read the probe QC CSV file into a data frame
failed.probes <- read.csv(opt$flagged, header = TRUE)[[1]]

### Filter Each R Data Object

# --- Function to load, filter, and save an RData object ---
process_minfi_object <- function(file_path, object_name, failed.probes, base_suffix) {
  if (!is.null(file_path) && file.exists(file_path)) {
    message(paste("Loading", object_name, "from:", file_path))
    # Load the R data object. It will create a variable with 'object_name' (e.g., 'MSet').
    load(file_path)

    # Get the loaded object dynamically by its name.
    current_object <- get(object_name)

    # Filter the object by removing failed probes
    initial_probe_count <- nrow(current_object)
    probes_to_remove_indices <- which(rownames(current_object) %in% failed.probes)
    current_object_filtered <- current_object[-probes_to_remove_indices,]
    
    message(paste("Filtered", object_name, ": Removed", 
                  initial_probe_count - nrow(current_object_filtered), 
                  "probes out of", initial_probe_count))

    # Subset the object to keep only the desired samples.
    # Assign it directly back to the original object_name variable.
    assign(object_name, current_object_filtered)

    # Construct the output file name.
    base_file_name <- tools::file_path_sans_ext(basename(file_path))
    output_file_name <- file.path(dirname(file_path), paste0(base_file_name, base_suffix, ".RData"))

    # Save the *modified* object (which has the original 'object_name' variable name).
    save(list = object_name, file = output_file_name)
    message(paste("  Filtered and saved:", output_file_name))

    message(paste("Filtered and saved:", output_file_name))
  } else if (!is.null(file_path)) {
    message(paste("Warning: File not found for", object_name, ":", file_path, ". Skipping filtering for this object."))
  } else {
    message(paste("No file path provided for", object_name, ". Skipping filtering for this object."))
  }
}


# --- Function to load, filter, and save an RData object ( RGset or RGsetEXT) ---
process_minfi_RGobject <- function(file_path, object_name, failed.probes, base_suffix) {
  if (!is.null(file_path) && file.exists(file_path)) {
    message(paste("Loading", object_name, "from:", file_path))
    # Load the R data object. It will create a variable with 'object_name' (e.g., 'MSet').
    load(file_path)

    # Get the loaded object dynamically by its name.
    current_object <- get(object_name)

    # Filter the object by removing failed probes
    initial_probe_count <- nrow(current_object)
    current_object_filtered = subsetByLoci(current_object, includeLoci = NULL, excludeLoci = failed.probes,
             keepControls = TRUE, keepSnps = TRUE)
    
    message(paste("Filtered", object_name, ": Removed", 
                  initial_probe_count - nrow(current_object_filtered), 
                  "probes out of", initial_probe_count))

    # Subset the object to keep only the desired samples.
    # Assign it directly back to the original object_name variable.
    assign(object_name, current_object_filtered)

    # Construct the output file name.
    base_file_name <- tools::file_path_sans_ext(basename(file_path))
    output_file_name <- file.path(dirname(file_path), paste0(base_file_name, base_suffix, ".RData"))

    # Save the *modified* object (which has the original 'object_name' variable name).
    save(list = object_name, file = output_file_name)
    message(paste("Filtered and saved:", output_file_name))

    message(paste("Filtered and saved:", output_file_name))
  } else if (!is.null(file_path)) {
    message(paste("Warning: File not found for", object_name, ":", file_path, ". Skipping filtering for this object."))
  } else {
    message(paste("No file path provided for", object_name, ". Skipping filtering for this object."))
  }
}


# --- Process each specified RData file ---

# Process MSet
process_minfi_object(opt$mset, "MSet", failed.probes, opt$base_suffix)

# Process grSet
process_minfi_object(opt$grset, "grSet", failed.probes, opt$base_suffix)

# For RGset and RGsetEXT, you cannot simply subset using subsetByLoci
# Process RGset
process_minfi_RGobject(opt$rgset, "RGset", failed.probes, opt$base_suffix)

# Process RGsetEXT
process_minfi_RGobject(opt$rgsetext, "RGsetEXT", failed.probes, opt$base_suffix)

message("\nSample filtering process completed for all specified R data files.")