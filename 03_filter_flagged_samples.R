# Author: Vartika Bisht
# Date: 24 May 2025
# Description: This script filters samples from multiple user-specified R data objects
# based on a list of flagged samples provided in a CSV file. Samples identified
# as "flagged" (where 'Total_Occurrences' is 2 or more) in the CSV are removed
# from each input R data object. The filtered objects are then saved with a
# user-defined suffix appended to their original file name, just before the '.Rdata' extension.
# Usage:
# Rscript 03_filter_flagged_samples.R -s sample_qc.csv -m MSet.Rdata -g grSet.Rdata -r ratioSet.Rdata -R RGset.Rdata -E RGsetEXT.Rdata -b _sample_filtered
# Input:
#   -s, --sample_qc_file: Path to the CSV file containing sample QC information (e.g., sample_qc.csv).
#                         This file is expected to have a column named 'Total_Occurrences'
#                         and a column named 'Sample' containing sample identifiers.
#   -m, --mset: Path to the MSet.Rdata file (optional, specify if you want to process this file).
#   -g, --grset: Path to the grSet.Rdata file (optional).
#   -r, --ratioset: Path to the ratioSet.Rdata file (optional).
#   -R, --rgset: Path to the RGset.Rdata file (optional).
#   -E, --rgsetext: Path to the RGsetEXT.Rdata file (optional).
#                  You can specify one or more of these R data file inputs.
#
# Output:
#   Filtered R data files, named by adding the specified base suffix before the '.Rdata' extension.
#   For example, if the input is `MSet.Rdata` and the base suffix is `_filtered`, the output file
#   will be `MSet_filtered.Rdata`. These files are saved in the same directory as their input counterparts.
#
# Parameters:
#   -b, --base_suffix: A character string to append to the object's original file name before the '.Rdata' extension
#                      when saving the filtered file (e.g., '_filtered', '_sample_cleaned').
#                      Defaults to '_sample_filtered'.
#   - The 'Total_Occurrences' threshold for flagging samples (>= 2) is hardcoded within the script.

library(optparse)
library(minfi) # Required for handling minfi objects like RGChannelSet, MethylSet, GenomicRatioSet

# Parse command line arguments
option_list <- list(
  make_option(c("-s", "--sample_qc_file"), type = "character", default = NULL,
              help = "Path to the sample QC CSV file (e.g., sample_qc.csv)"),
  make_option(c("-m", "--mset"), type = "character", default = NULL,
              help = "Path to the MSet.Rdata file"),
  make_option(c("-g", "--grset"), type = "character", default = NULL,
              help = "Path to the grSet.Rdata file"),
  make_option(c("-r", "--ratioset"), type = "character", default = NULL,
              help = "Path to the ratioSet.Rdata file"),
  make_option(c("-R", "--rgset"), type = "character", default = NULL,
              help = "Path to the RGset.Rdata file"),
  make_option(c("-E", "--rgsetext"), type = "character", default = NULL,
              help = "Path to the RGsetEXT.Rdata file"),
  make_option(c("-b", "--base_suffix"), type = "character", default = "_sample_filtered",
              help = "Suffix to add to output file names (e.g., '_filtered'). Default is '_sample_filtered'.")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Validate input argument: check if sample QC file path is provided and exists
if (is.null(opt$sample_qc_file) || !file.exists(opt$sample_qc_file)) {
  stop("Please provide the path to the sample QC CSV file using -s or --sample_qc_file option and ensure it exists.")
}

# Collect all R data file paths into a named list for easier iteration
# The name in the list should match the object name loaded from the .Rdata file
input_rdata_paths <- list(
  MSet = opt$mset,
  grSet = opt$grset,
  ratioSet = opt$ratioset,
  RGset = opt$rgset,
  RGsetEXT = opt$rgsetext
)

# Filter out NULL entries for any R data files not specified by the user
input_rdata_paths <- input_rdata_paths[!sapply(input_rdata_paths, is.null)]

if (length(input_rdata_paths) == 0) {
  stop("No R data files were provided. Please specify at least one R data file using -m, -g, -r, -R, or -E options.")
}

# Check if all specified R data files actually exist on disk
for (file_path in input_rdata_paths) {
  if (!file.exists(file_path)) {
    stop(paste("The specified R data file does not exist:", file_path))
  }
}


### Identify Samples to Exclude

# Read the sample QC CSV file into a data frame
flagged_samples_df <- read.csv(opt$sample_qc_file, header = TRUE, stringsAsFactors = FALSE)

# Filter the data frame to get samples where 'Total_Occurrences' is 2 or more
flagged_samples_to_exclude <- flagged_samples_df[flagged_samples_df$Total_Occurrences >= 2, ]$Sample

# Provide a summary message about the samples identified for exclusion
if (length(flagged_samples_to_exclude) == 0) {
  message("No samples were flagged for exclusion based on 'Total_Occurrences' >= 2.")
} 

### Filter Each R Data Object

# Iterate through each specified R data file
for (object_name in names(input_rdata_paths)) {
    file_path <- input_rdata_paths[[object_name]]
    message(paste("\nProcessing file:", file_path))

    # Load the R data object into the current R environment.
    # This will create a variable (e.g., 'MSet', 'grSet') with the name 'object_name' in your R session.
    load(file_path)

    # Get the loaded object dynamically by its name to perform operations
    current_object <- get(object_name)

    # Determine the indices of samples to keep (i.e., those not present in the exclusion list)
    samples_to_keep_idx <- which(!current_sample_names %in% flagged_samples_to_exclude)

    # Subset the object to keep only the desired samples
    filtered_object <- current_object[, samples_to_keep_idx]

    # Assign the filtered object back to its original name in the environment.
    # This is crucial so that `save(list = object_name, ...)` correctly saves the modified object.
    assign(object_name, filtered_object)

    # Construct the output file name by inserting the base suffix just before the '.Rdata' extension.
    # `tools::file_path_sans_ext` gets the filename without extension (e.g., "MSet").
    base_file_name <- tools::file_path_sans_ext(basename(file_path))
    output_file_name <- file.path(dirname(file_path), paste0(base_file_name, opt$base_suffix, ".Rdata"))

    # Save the filtered object to the new R data file
    save(list = object_name, file = output_file_name)
    message(paste("  Filtered and saved:", output_file_name))

    # Remove the loaded object from the environment tosave space
    # when loading the next file in the loop.
    rm(list = object_name)
    gc() # Clean up memory
}

message("\nSample filtering process completed for all specified R data files.")