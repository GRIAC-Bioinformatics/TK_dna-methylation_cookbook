#!/usr/bin/env Rscript
# Date: 30 May 2020
# Author: Vartika Bisht
#
# Intensity Comparison for Methylation Array Quality Control
#
# This script is designed to help you visually assess the quality of samples
# from your methylation array dataset. It does this by plotting the
# log median methylated intensity against the log median unmethylated intensity.
# The QC information for this plot is extracted from an MSet object using
# the minfi::getQC() function. This function provides a DataFrame
# with two crucial columns: mMed (the chip-wide median of the methylated
# channel intensities) and uMed (the chip-wide median of the unmethylated
# channel intensities). Here's how the dataset looks like:
#
# DataFrame with 6 rows and 2 columns
#                         mMed      uMed
#                       <numeric> <numeric>
# 6264488065_R01C01   11.64746   11.15292
# 6264488065_R02C01    9.28771    9.09803
# 6264488065_R03C01   11.52454   11.12993
# 6264488065_R04C01   11.93995   11.28598
# 6264488065_R05C01   11.98655   11.31288
# 6264488065_R06C01   12.00633   11.35094
#
# Important Note: The cutoff threshold used to identify "bad" samples is
# an empirical value. It's arbitrary and might need to be adjusted based on
# the specific characteristics of your dataset to accurately reflect sample quality.
#
# ---
#
# ### Command Line Parameters
#
# This script is designed to be run from the command line and accepts the
# following arguments:
#
# * -m, --mset (Required): The file path to your MSet.RData file.
# * -c, --cutoff (Optional): A numeric value that sets the intensity
#   cutoff. Samples where both the methylated and unmethylated median
#   intensities fall below this number will be flagged. The default value is 10.5.
# * -f, --flagged (Required): The file path for the CSV output.
#   This file will list the names of any samples that have been flagged as low quality.
# * -p, --pdf (Required): The file path for the output PDF. This is
#   where the generated QC plot will be saved
#
# ---
#
# ### Usage Example : Rscript 02_intensity_flagged_samples.R --mset MSet.RData --cutoff 10.5 --flagged Flagged_intensity.csv --pdf Intensity_flagged_samples.pdf

# Silent library loading and error handling
suppressPackageStartupMessages({
  tryCatch({
  library(minfi)
  library(optparse)
  }, error = function(e) {
    stop("Package loading failed: ", e$message, call. = FALSE)
  })
})


# Set up command line options
option_list <- list(
  make_option(c("-m", "--mset"), type="character", default=NULL,
              help="Path to MSet.RData file", metavar="FILE"),
  make_option(c("-c", "--cutoff"), type="numeric", default=10.5,
              help="Cutoff for methylated and unmethylted intensities", metavar="NUM"),
  make_option(c("-f", "--flagged"), type="character", default=NULL,
              help="CSV file to save samples which have been flagged", metavar="FILE"),
  make_option(c("-p", "--pdf"), type="character", default=NULL,
              help="Output PDF file path", metavar="FILE")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$mset) || !file.exists(opt$mset)) {
  stop("Please provide a valid path to the MSet.RData file using -m or --mset option.")
}
if (is.null(opt$flagged)) {
  stop("Please provide the path to save flagged samples using -f or --flagged option.")
}
if (is.null(opt$pdf)) {
  stop("Please provide the path to save the PDF output using -p or --pdf option.")
}

# Load the MSet object
load(opt$mset)

message("Plotting log median (met) vs log median(un-met)...")

tryCatch({
  qc <- minfi::getQC(MSet)

  # Calculate values
  meds <- (qc$mMed + qc$uMed) / 2
  bad_sample_idx <- intersect(which(qc$mMed < opt$c), which(qc$uMed < opt$c))

  write.csv(data.frame(intersect(names(which(qc$mMed < 10.5)), names(which(qc$uMed < 10.5)))),
    file = opt$flagged, row.names = FALSE
  )

  pdf(
    file = opt$pdf,
    width = 8, height = 6
  )

  # Set up plot with better margins
  par(mar = c(5, 5, 4, 2) + 0.1)

  # Main plot with enhanced title
  plot(qc$mMed, qc$uMed,
    xlab = "Methylated median intensity (log2)",
    ylab = "Unmethylated median intensity (log2)",
    col = ifelse(1:nrow(qc) %in% bad_sample_idx, "red", "black"),
    pch = ifelse(1:nrow(qc) %in% bad_sample_idx, 17, 16),
    cex = 1.2,
    xlim = c(5, 14),
    ylim = c(5, 14),
    main = "Per Sample Median Intensities",
    cex.main = 1.1,
    las = 1
  )

  # Add DEFAULT cutoff lines (blue)
  abline(h = opt$c, lty = 3, col = "blue")
  abline(v = opt$c, lty = 3, col = "blue")

  # Improved legend
  legend("topleft",
    legend = c(
      "Good samples",
      "Flagged samples",
      sprintf("Current cutoff (%.1f)", opt$cutoff)
    ),
    pch = c(16, 15, NA),
    col = c("black", "red", "blue"),
    lty = c(NA, NA, 3),
    bty = "n",
    cex = 0.8
  )

  # Add grid
  grid(col = "gray90", lty = "dotted")
  dev.off()
  message("QC plot and flagged samples file generated successfully.")

}, error = function(e) {
  message("An error occurred: ", e$message)
  # You might want to add more specific error handling or logging here
})

