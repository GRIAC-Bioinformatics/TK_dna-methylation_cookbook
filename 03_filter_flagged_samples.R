# Author: Vartika Bisht
# Date: 24 May 2025
# Description: This script filters samples from multiple user-specified R data objects
# based on a list of flagged samples provided in a CSV file. Samples identified
# as "flagged" (where 'Total_Occurrences' is 2 or more) in the CSV are removed
# from each input R data object. The filtered objects are then saved with a
# user-defined suffix appended to their original file name, just before the '.RData' extension.
# Usage:
# Rscript 03_filter_flagged_samples.R -s sample_qc.csv -m MSet.RData -g grSet.RData -r ratioSet.RData -R RGset.RData -E RGsetEXT.RData -b _sample_filtered
# Input:
#   -s, --sample_qc_file: Path to the CSV file containing sample QC information (e.g., sample_qc.csv).
#                         This file is expected to have a column named 'Total_Occurrences'
#                         and a column named 'Sample' containing sample identifiers.
#   -m, --mset: Path to the MSet.RData file (optional, specify if you want to process this file).
#   -g, --grset: Path to the grSet.RData file (optional).
#   -r, --ratioset: Path to the ratioSet.RData file (optional).
#   -R, --rgset: Path to the RGset.RData file (optional).
#   -E, --rgsetext: Path to the RGsetEXT.RData file (optional).
#                  You can specify one or more of these R data file inputs.
#
# Output:
#   Filtered R data files, named by adding the specified base suffix before the '.RData' extension.
#   For example, if the input is `MSet.RData` and the base suffix is `_filtered`, the output file
#   will be `MSet_filtered.RData`. These files are saved in the same directory as their input counterparts.
#
# Parameters:
#   -b, --base_suffix: A character string to append to the object's original file name before the '.RData' extension
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
              help = "Path to the MSet.RData file"),
  make_option(c("-g", "--grset"), type = "character", default = NULL,
              help = "Path to the grSet.RData file"),
  make_option(c("-r", "--ratioset"), type = "character", default = NULL,
              help = "Path to the ratioSet.RData file"),
  make_option(c("-R", "--rgset"), type = "character", default = NULL,
              help = "Path to the RGset.RData file"),
  make_option(c("-E", "--rgsetext"), type = "character", default = NULL,
              help = "Path to the RGsetEXT.RData file"),
  make_option(c("-b", "--base_suffix"), type = "character", default = "_sample_filtered",
              help = "Suffix to add to output file names (e.g., '_filtered'). Default is '_sample_filtered'."),
  make_option(c("-p", "--output_plot_name"), type = "character", default = NULL,
              help = "Name for the output PDF plot file containing all plots (required)."),
  make_option(c("-f", "--min_flag_overlap"), type = "character", default = NULL,
              help = "Minimum number of overlapping flags required to filter a sample (required).")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Validate input argument: check if sample QC file path is provided and exists
if (is.null(opt$sample_qc_file) || !file.exists(opt$sample_qc_file)) {
  stop("Please provide the path to the sample QC CSV file using -s or --sample_qc_file option and ensure it exists.")
}

# Collect all R data file paths into a named list for easier iteration
# The name in the list should match the object name loaded from the .RData file
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


# ### Identify Samples to Exclude

# # Read the sample QC CSV file into a data frame
# flagged_samples_df <- read.csv(opt$sample_qc_file, header = TRUE, stringsAsFactors = FALSE)

# # Filter the data frame to get samples where 'Total_Occurrences' is 2 or more
# flagged_samples_to_exclude <- flagged_samples_df[flagged_samples_df$Total_Occurrences >= as.numeric(opt$min_flag_overlap), ]$Sample

# # Provide a summary message about the samples identified for exclusion
# if (length(flagged_samples_to_exclude) == 0) {
#   message("No samples were flagged for exclusion based on 'Total_Occurrences' >= ",min_flag_overlap,".")
# } 

# ### Filter Each R Data Object

# # Iterate through each specified R data file
# for (object_name in names(input_rdata_paths)) {
#     file_path <- input_rdata_paths[[object_name]]
#     message(paste("\nProcessing file:", file_path))

#     # Load the R data object into the current R environment.
#     # This will create a variable (e.g., 'MSet', 'grSet') with the name 'object_name' in your R session.
#     load(file_path)

#     # Get the loaded object dynamically by its name to perform operations
#     current_object <- get(object_name)

#     # Determine the indices of samples to keep (i.e., those not present in the exclusion list)
#     samples_to_keep_idx <- which(!sampleNames(current_object) %in% flagged_samples_to_exclude)

#     # Subset the object to keep only the desired samples
#     filtered_object <- current_object[, samples_to_keep_idx]

#     # Construct the output file name by inserting the base suffix just before the '.RData' extension.
#     # `tools::file_path_sans_ext` gets the filename without extension (e.g., "MSet").
#     base_file_name <- tools::file_path_sans_ext(basename(file_path))
#     output_file_name <- file.path(dirname(file_path), paste0(base_file_name, opt$base_suffix, ".RData"))

#     # Save the filtered object to the new R data file
#     save(filtered_object, file = output_file_name)
#     message(paste("  Filtered and saved:", output_file_name))

#     # Remove the loaded object from the environment tosave space
#     # when loading the next file in the loop.
#     rm(filtered_object)
#     rm(list = object_name)
#     gc() # Clean up memory
# }

# message("\nSample filtering process completed for all specified R data files.")

## Overall distributions of Beta values for each sample
# Here we expect beta values presenting values close to zero or one. 
# This is an indication of unmethylatilated CpG sites when close to zero and methylated when close to one.
# if you have NAs due to 0 intesity values you will have errors here so to avoid this, use offset

tryCatch({
    message("Plotting overall distributions of Beta values for each sample..")

    # Open PDF device for plotting
    pdf(opt$output_plot_name, width = 8, height = 6)
    message(paste("PDF plot will be saved to:", opt$output_plot_name))

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
    filtered_mset_path <- file.path(dirname(file_path), paste0(tools::file_path_sans_ext(basename(opt$mset)), opt$base_suffix, ".RData"))
    message(paste("Loading filtered MSet object from:", filtered_mset_path))

    load(filtered_mset_path) 

    # Will be loaded as filtered_object 
    if (!exists("filtered_object")) {
        stop("Error: 'filtered_object' not found after loading filtered object file.")
    }

    # Rename filtered_object to MSet
    MSet_filtered <- filtered_object

    # Remove the intermediate filtered_object to free memory
    rm(filtered_object)

    # Extract beta values for the filtered MSet
    beta_filtered_matrix <- getBeta(MSet_filtered, offset = 100)
    message("Successfully extracted beta values from filtered MSet.")

    # Remove the filtered MSet object as its beta values are now extracted
    rm(MSet_filtered)

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