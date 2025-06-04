#!/usr/bin/env Rscript
# Date: 30 May 2020
# Author: Vartika Bisht
#
# ## Detection P-value Filtering
#
# This script performs quality control by assessing detection p-values, which are crucial for identifying unreliable samples in methylation array data.
#
# The minfi::detectionP() function calculates these p-values by comparing the total signal
# (methylated + unmethylated) for each probe against the background signal level, which is estimated
# from negative control probes defined by Illumina. Small p-values indicate a reliable signal, while
# large p-values (typically greater than 0.01) suggest a poor-quality or unreliable signal. We 
# visualise the samples based on their mean detection p-values across all probes and flag any samples
# that exceed a specified cutoff.
#
# Important Note: The choice of p-value cutoff is arbitrary and should be determined based on your
# specific dataset and research goals.
#
# ---
#
# ### Command Line Parameters
#
# This script accepts the following command-line arguments:
#
# * -r, --rgset (Required): The file path to your RGset object. 
# * -f, --flagged (Required): The file path for the CSV output. This file will list the names
#   of any probes that have been flagged due to high detection p-values.
# * -c, --cutoff (Optional): A numeric value that sets the detection p-value cutoff. Probes
#   with p-values greater than this number will be considered unreliable. The default value is 0.01.
# * -p, --pdf (Required): The file path for the output PDF.
#
# ### Usage : Rscript 02_sample_detection_p_value.R --rgset RGset.RData --flagged Flagged_detectionP.csv --cutoff 0.01 --pdf detectionP_flagged_samples.pdf
#

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
    make_option(c("-r", "--rgset"), type="character", default=NULL,
                help="Path to RGset file", metavar="FILE"),
    make_option(c("-f", "--flagged"), type="character", default=NULL,
                help="CSV file to save flagged probes which have high detection p-values", metavar="FILE"),
    make_option(c("-c", "--cutoff"), type="numeric", default=0.01,
                help="Cutoff for detection p-value (default: 0.01)", metavar="NUM"),
    make_option(c("-p", "--pdf"), type="character", default=NULL,
                help="Output PDF file path", metavar="FILE")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$rgset) || !file.exists(opt$rgset)) {
  stop("Please provide a valid path to the RGset file using -r or --rgset option.")
}
if (is.null(opt$flagged)) {
  stop("Please provide the path to save flagged probes using -f or --flagged option.")
}
if (is.null(opt$pdf)) {
  stop("Please provide the path to save the PDF output using -p or --pdf option.")
}
tryCatch({
  # Load the RGset object
  message("Loading RGset object...")
  load(opt$rgset) # Ensure opt$rgset is defined, e.g., from optparse

  message("Calculating detection p-values...")
  # An error occurred: Argument 'useNames' must be either TRUE or FALSE
  # detP <- minfi::detectionP(RGset)
  # That is why we need to the source code
  # Extract data to pass to low-level function that constructs `detP`
  locusNames <- getManifestInfo(RGset, "locusNames")
  controlIdx <- getControlAddress(RGset, controlType = "NEGATIVE")
  Red <- getRed(RGset)
  Green <- getGreen(RGset)
  TypeI.Red <- getProbeInfo(RGset, type = "I-Red")
  TypeI.Green <- getProbeInfo(RGset, type = "I-Green")
  TypeII <- getProbeInfo(RGset, type = "II")
  # Set up output matrix with appropriate dimensions and type
  detP <- matrix(NA_real_,
                  nrow = length(locusNames),
                  ncol = ncol(Red),
                  dimnames = list(locusNames, colnames(Red)))

  # Compute summary statistics needed for calculations
  rBg <- Red[controlIdx, , drop = FALSE]
  rMu <- colMedians(rBg, na.rm = TRUE, useNames = TRUE)
  rSd <- matrixStats::colMads(rBg)
  gBg <- Green[controlIdx, , drop = FALSE]
  gMu <- colMedians(gBg , na.rm = TRUE , useNames = TRUE)
  gSd <- matrixStats::colMads(gBg)

  # Fill output matrix
  for (j in seq_len(ncol(detP))) {
      # Type I Red
      intensity <- Red[TypeI.Red$AddressA, j] + Red[TypeI.Red$AddressB, j]
      detP[TypeI.Red$Name, j] <- pnorm(
          q = intensity,
          mean = 2 * rMu[j],
          sd = 2 * rSd[j],
          lower.tail = FALSE)
      # Type I Green
      intensity <- Green[TypeI.Green$AddressA, j] +
          Green[TypeI.Green$AddressB, j]
      detP[TypeI.Green$Name, j] <- pnorm(
          q = intensity,
          mean = 2 * gMu[j],
          sd = 2 * gSd[j],
          lower.tail = FALSE)
      # Type II
      intensity <- Red[TypeII$AddressA, j] + Green[TypeII$AddressA, j]
      detP[TypeII$Name, j] <- pnorm(
          q = intensity,
          mean = rMu[j] + gMu[j],
          sd = rSd[j] + gSd[j],
          lower.tail = FALSE)
  }

  message("Detection p-values calculated successfully.")

  # Calculate mean detection p-values for each sample
  mean_detP <- colMeans(detP, na.rm = TRUE)

  # Convert mean detection p-values to -log10 scale
  neglog10_mean_detP <- -log10(mean_detP)

  # Define the -log10(opt$cutoff) cutoff
  cutoff <- -log10(as.numeric(opt$cutoff))

  # Flag samples with mean detection p-value > opt$cutoff
  flagged_sample_names <- names(which(mean_detP > as.numeric(opt$cutoff)))
  num_flagged_samples <- length(flagged_sample_names)

  # Create a data frame for flagged samples and save to CSV
  if (num_flagged_samples > 0) {
    flagged_sample_df <- data.frame(Sample = flagged_sample_names)
    write.csv(flagged_sample_df, file = opt$flagged, row.names = FALSE)
    message(paste0("Flagged ", num_flagged_samples, " samples with mean p-value > ",opt$cutoff," saved to ", opt$flagged))
  } else {
    message("No samples flagged with mean p-value > 0.01.")
    # Ensure an empty file is still created if no samples are flagged
    write.csv(data.frame(Sample = character()), file = opt$flagged, row.names = FALSE)
  }

  # Open PDF for plotting
  pdf(opt$pdf, width = 10, height = 7) # Increased width for sample names on scatter plot

  # --- PLOT 1: Scatter plot of mean detection p-values ---
  message("Plotting mean detection p-values (scatter plot)...")

  # Set up plot margins to accommodate sample names
  par(mar = c(8, 4, 4, 2) + 0.1) # Bottom margin increased for x-axis labels

  plot(
    x = 1:length(neglog10_mean_detP),
    y = neglog10_mean_detP,
    type = "n", # Do not draw points initially
    xlab = "", # No generic X-axis label, will use sample names
    ylab = "-log10 Mean Detection P-value",
    main = paste0("Mean Detection P-value Per Sample (Red line: p = ",opt$cutoff,")"),
    xaxt = "n", # Suppress default x-axis
  )

  # Add points for samples not flagged
  points(
    x = which(!(names(neglog10_mean_detP) %in% flagged_sample_names)),
    y = neglog10_mean_detP[!(names(neglog10_mean_detP) %in% flagged_sample_names)],
    pch = 16, # Solid circle
    col = "black",
    cex = 0.8
  )

  # Add sample names for flagged samples
  if (num_flagged_samples > 0) {
      points(
    x = which((names(neglog10_mean_detP) %in% flagged_sample_names)),
    y = neglog10_mean_detP[(names(neglog10_mean_detP) %in% flagged_sample_names)],
    pch = 16, # Solid circle
    col = "red",
    cex = 0.8
  )
  }

  # Add red line at -log10(opt$cutoff)
  abline(h = cutoff, col = "red", lty = 2, lwd = 1)

  # Add custom x-axis labels (sample names)
  axis(
    side = 1,
    at = 1:length(neglog10_mean_detP),
    labels = names(neglog10_mean_detP),
    las = 2, # Rotate labels vertically
    cex.axis = 0.7 # Adjust label size
  )

  dev.off()
  message("Detection p-value plots generated successfully and saved to ", opt$pdf)

}, error = function(e) {
  message("An error occurred: ", e$message)
  # Ensure dev.off() is called if pdf() was opened to avoid corrupted files
  if (!is.null(dev.list())) {
    dev.off()
  }
})