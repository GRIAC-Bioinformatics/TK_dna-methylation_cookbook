#!/usr/bin/env Rscript
# =============================================================================
# 03_intensity_flagged_samples.R
# Author: Vartika Bisht
# Created: 13-Aug-2025
# Description:
# This script performs a quality control check on DNA methylation array data by
# plotting the log median methylated intensity against the log median
# unmethylated intensity. Samples falling below a user-defined threshold are
# identified and flagged as potential outliers due to low signal intensity.
# Arguments:
#    --mset                 → Path to the MSet object file (.RData)
#    --cutoff               → Numeric cutoff for methylated and unmethylated intensities (default: 10.5)
#    --flagged              → Path to a CSV file to save the list of flagged sample IDs
#    --pdf_output           → Output PDF file path for the intensity plot
# Usage:
#   Rscript 03_intensity_flagged_samples.R --mset <FILE.RData> \
#       --cutoff 10.5 --pdf_output <output.pdf> --flagged <flagged_samples.csv>
# Notes:
#   - The cutoff value is empirical and may need to be adjusted based on the
#     specific characteristics of your dataset.
#   - Samples below the cutoff on either axis are considered low quality and
#     are recommended for removal from downstream analysis.
# =============================================================================

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
  make_option(c("-p", "--pdf_output"), type="character", default=NULL,
              help="Output PDF file path", metavar="FILE")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)


# VALIDATE REQUIRED ARGUMENTS
required_args <- c("mset", "cutoff", "flagged","pdf_output")

for (arg in required_args) {
  if (is.null(opt[[arg]])) {
    stop(paste("Error: Required argument --", arg, " is missing. Use --help for more information.", sep = ""))
  }
}


# Load the MSet object
load(opt$mset)

message("Plotting log median (met) vs log median(un-met)...")

tryCatch({
  message("Generating QC dataframe using colMedians of methylated and umethylated signals..")
  message("Summarise all probe signals from one sample to median..")
  qc <- DataFrame(mMed = log2(colMedians(getMeth(MSet), na.rm = TRUE, useNames = TRUE)),
                  uMed = log2(colMedians(getUnmeth(MSet), na.rm = TRUE, useNames = TRUE)))

  rownames(qc) <- colnames(MSet)

  # Calculate values
  message("Intensity cutoff of ",opt$cutoff," provided.")
  bad_sample_idx <- intersect(which(qc$mMed < opt$cutoff), which(qc$uMed < opt$cutoff))

  message(length(bad_sample_idx) ," samples flagged due to their median intensity being < ",opt$cutoff," in both methylated and unmethylated channel.")
  message("Writing these sample IDs down in a CSV file : ",opt$flagged)
  write.csv(data.frame(fagged = intersect(names(which(qc$mMed < opt$cutoff)), names(which(qc$uMed < opt$cutoff)))),
    file = opt$flagged, row.names = FALSE
  )

  pdf(
    file = opt$pdf_output,
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
  abline(h = opt$cutoff, lty = 3, col = "blue")
  abline(v = opt$cutoff, lty = 3, col = "blue")

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

