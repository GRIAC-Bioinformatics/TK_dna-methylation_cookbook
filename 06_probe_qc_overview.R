#!/usr/bin/env Rscript
# =============================================================================
# 06_probe_qc_overview.R
# Author: Vartika Bisht
# Created: 13-Aug-2025
# Description:
# This script generates a comprehensive overview of all probes flagged during
# various upstream quality control steps. It combines flagging information from
# multiple QC analyses into a single summary report, providing plots, such as an
# upset plot, and a combined data file to identify probes with significant
# quality issues.
# Arguments:
#    --flaggedlist           → CSV file with a list of QC flags and the paths to their respective flagged probe files
#    --min_flag_overlap      → Minimum number of flags a probe must have to be included in the upset plot and summary
#    --flaggedcombined       → Output CSV file with all combined flag information
#    --pdf_output            → Output PDF file containing all QC overview plots
# Usage:
#   Rscript 06_probe_qc_overview.R --flaggedlist <flag_paths.csv> \
#     --min_flag_overlap 2 --flaggedcombined <combined_flags.csv> \
#     --pdf_output <probe_qc_overview.pdf>
# Notes:
#   - This script integrates flags from various quality control scripts that target probes,
#     including like , but not restricted to  :
#      - High Detection P-values
#      - Low Bead Counts
#      - SNP-containing Probes
#      - High Intensity Probes
#      - Sex Chromosome Location
#      - Cross-reactivity
#   - The combined flag information is used to generate an upset plot , which visualizes the
#     intersection of different probe flag categories. This helps identify probes that fail
#     multiple QC metrics and assess the overlap between different filtering criteria.
#   - The min_flag_overlap parameter allows for focusing the analysis on probes with more
#     severe, overlapping quality issues.
# =============================================================================

suppressPackageStartupMessages({
  tryCatch({
library(optparse)
library(ggplot2)
library(minfi)
library(UpSetR)
  }, error = function(e) {
    stop("Package loading failed: ", e$message, call. = FALSE)
  })
})

option_list <- list(
  make_option(c("-f", "--flaggedlist"), type = "character", default = NULL,
              help = "CSV file with list of QC flag and the path to the CSV file with the list of probes that failed that QC step (required)."),
  make_option(c("-m", "--min_flag_overlap"), type = "character", default = NULL,
              help = "Minimum number of overlapping flags required to filter a probes (required)."),
  make_option(c("-o", "--flaggedcombined"), type = "character", default = NULL,
              help = "Name for the output CSV file (required)."),
  make_option(c("-p", "--pdf_output"), type = "character", default = NULL,
              help = "Name for the output PDF plot file containing all plots (required).")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)


# VALIDATE REQUIRED ARGUMENTS
required_args <- c( "flaggedcombined" , "flaggedlist" , "pdf_output", "min_flag_overlap")

for (arg in required_args) {
  if (is.null(opt[[arg]])) {
    stop(paste("Error: Required argument --", arg, " is missing. Use --help for more information.", sep = ""))
  }
}


message("Load CSV file with the list of QC flags and the probes that failed each QC step from ", opt$flaggedlist, " ...")
flaggedlist <- tryCatch({
  read.csv(opt$flaggedlist, header = FALSE)
}, error = function(e) {
  stop("Failed to read flagged list file: ", e$message)
})

# The first column of flagged.list is the QC flag name, and the second column is the path to the CSV file with the list of probes that failed that QC step
# Make a list of list with the list name as the QC flag name and the list of probes that failed that QC step
flagged <- list()
for (i in 1:nrow(flaggedlist)) {
  flag <- flaggedlist$V1[i]
  path <- flaggedlist$V2[i]
  message("Reading flagged probes for flag: ", flag)
  flag_probes <- tryCatch({
    read.csv(path, header = TRUE)[[1]]
  }, error = function(e) {
    stop("Failed to read flagged probes file for flag '", flag, "': ", e$message)
  })

  if(length(flag_probes) == 0) {
    message("No probes found for flag: ", flag)
  } else {
    message("Found ", length(flag_probes), " probes for flag: ", flag)
    flagged[[flag]] <- flag_probes
  }
}

message("Prosessing flags : ", paste(names(flagged), collapse = " , "))


# Create a binary data.frame with probes as rows and flags as columns
# Each cell will be TRUE if the probe belongs to that flag, FALSE otherwise
all_flagged_probes <- unique(unlist(flagged))
all_flags <- as.character(names(flagged))

# Initialize empty data.frame of FALSE
flagged_df <- as.data.frame(
  matrix(0, nrow = length(all_flagged_probes), ncol = length(all_flags)),
  row.names = all_flagged_probes
)
colnames(flagged_df) <- all_flags

# Fill in TRUE where probe belongs to a flag
for (flag in all_flags) {
  message("Reading flagged probes for flag: ", flag)
  probelist <- as.character(flagged[[flag]])
  flagged_df[probelist, flag] <- 1
}

# Add 'Total Occurrences' column (based on combined categories)
flagged_df$Total_Occurrences <- rowSums(flagged_df)
flagged_df$Probes <- rownames(flagged_df)

# Sort the table by 'Total_Occurrences' in descending order
flagged_df <- flagged_df[order(flagged_df$Total_Occurrences, decreasing = TRUE),c("Probes", "Total_Occurrences", all_flags)]

write.csv(flagged_df, file = opt$flaggedcombined, row.names = FALSE) 
message(paste("Probe QC overview saved to:", opt$flaggedcombined))



# Plotting Section 
pdf(opt$pdf_output, width = 12, height = 8) # Open PDF device for all plots

# Apply `as.numeric` to convert these columns from character/logical to numeric (0 or 1)
flagged_df[all_flags] <- lapply(flagged_df[all_flags], as.numeric)

message("Generating UpSet plot..")

# --- Generate the UpSet plot ---
# The nsets parameter should now be dynamic based on the filtered list
tryCatch({
  if (nrow(flagged_df) > 0 && length(all_flags) > 0) {
    # Specify which columns in your data frame are the 'sets' for the UpSet plot
    sets_to_plot <- all_flags

    # Create the UpSet plot
    print(upset(
      flagged_df,               # Your data frame containing 0s and 1s for flags
      sets = sets_to_plot,        # The names of the columns to treat as sets
      main.bar.color = "steelblue", # Color for the main intersection bars (top)
      sets.bar.color = "darkgreen", # Color for the individual set size bars (left)
      matrix.color = "darkblue",    # Color for the dots in the intersection matrix
      point.size = 2.5,             # Size of the dots in the matrix
      line.size = 1,                # Thickness of the lines connecting dots in the matrix
      text.scale = c(1.3, 1.3, 1, 1, 1.5, 1.2), # Adjust font sizes for different plot elements
      order.by = "freq",            # Order intersections by their frequency (most common first)
      empty.intersections = "on",   # Include combinations that have zero probes
      mainbar.y.label = "Intersection Size (Number of Probes)", # Y-axis label for top bars
      sets.x.label = "Set Size (Total Probes per Flag)" # X-axis label for left bars
    ))

    message("UpSet Plot generated successfully.")

  } else {
    warning("No probes or flag columns found for UpSet plot.")
    message("Please ensure 'flagged_df' has rows and relevant flag columns (e.g., 'Flagged_...csv').")
  }
}, error = function(e) {
  message("Error generating UpSet Plot: ", e$message)
  message("Please check your data and ensure 'UpSetR' is installed and loaded correctly.")
})

dev.off()