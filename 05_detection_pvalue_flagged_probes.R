#!/usr/bin/env Rscript
# =============================================================================
# 05_detection_pvalue_flagged_probes.R
# Author: Vartika Bisht
# Created: 13-Aug-2025
# Description:
# This script performs quality control by filtering out unreliable probes based on
# their detection p-values. It identifies probes with high p-values across a
# significant fraction of samples, which may indicate systematic issues such as
# low intensity or scanner artifacts.
# Arguments:
#    --rgset                → Path to the RGChannelSet object file (.RData)
#    --flagged              → Output CSV file to list flagged probes
#    --cutoff               → Numeric cutoff for detection p-value (default: 0.01)
#    --threshold            → Fraction of samples a probe must pass the p-value cutoff in to be retained (default: 0.99)
#    --pdf_output           → Output PDF file path for QC plots
# Usage:
#   Rscript 05_detection_pvalue_flagged_probes.R --rgset <FILE.RData> \
#     --flagged <flagged_probes.csv> --cutoff 0.01 --threshold 0.99 \
#     --pdf_output <output.pdf>
# Notes:
#   - This script filters probes, not samples, based on detection p-values.
#   - A probe is flagged and removed if its detection p-value is greater than the
#     specified cutoff in more than 1 - threshold fraction of the samples.
#     For example, with a cutoff of 0.01 and a threshold of 0.99, probes with
#     p-values > 0.01 in more than 1% of the samples will be removed.
#   - Filtering unreliable probes is a critical step to ensure that downstream
#     analysis is based on robust and confident methylation measurements.
# =============================================================================

# Silent library loading and error handling
suppressPackageStartupMessages({
  tryCatch({
  library(minfi)
  library(optparse)
  library(ggplot2)
  library(matrixStats)
  library(reshape2)
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
    make_option(c("-t", "--threshold"), type="numeric", default=0.99,
                help="Fraction of samples in which the probe should pass the p value cut off, for a probe to pass detectionP filter", metavar="NUM"),          
    make_option(c("-p", "--pdf_output"), type="character", default=NULL,
                help="Output PDF file path", metavar="FILE")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

# VALIDATE REQUIRED ARGUMENTS
required_args <- c("rgset", "cutoff", "flagged" , "threshold" , "pdf_output")

for (arg in required_args) {
  if (is.null(opt[[arg]])) {
    stop(paste("Error: Required argument --", arg, " is missing. Use --help for more information.", sep = ""))
  }
}

load(opt$rgset)
# This RGset object has been filtered for "bad samples" in the previous step
# RGset --> this is expected to be the name of the loaded object
message("Loading sample RGset object ", opt$rgset, "...")

tryCatch({

  message("Starting detection p-value calculation...")

  # Extract data to pass to low-level function that constructs `detP`
  message("Extracting manifest and probe information...")
  locusNames <- getManifestInfo(RGset, "locusNames")
  controlIdx <- getControlAddress(RGset, controlType = "NEGATIVE")
  Red <- getRed(RGset)
  Green <- getGreen(RGset)
  TypeI.Red <- getProbeInfo(RGset, type = "I-Red")
  TypeI.Green <- getProbeInfo(RGset, type = "I-Green")
  TypeII <- getProbeInfo(RGset, type = "II")

  # Set up output matrix
  message("Setting up output matrix...")
  detP <- matrix(NA_real_,
    nrow = length(locusNames),
    ncol = ncol(Red),
    dimnames = list(locusNames, colnames(Red))
  )

  # Compute summary statistics
  message("Computing summary statistics...")
  rBg <- Red[controlIdx, , drop = FALSE]
  rMu <- colMedians(rBg, na.rm = TRUE, useNames = TRUE)
  rSd <- matrixStats::colMads(rBg)

  gBg <- Green[controlIdx, , drop = FALSE]
  gMu <- colMedians(gBg, na.rm = TRUE, useNames = TRUE)
  gSd <- matrixStats::colMads(gBg)

  # Fill output matrix
  message("Filling output matrix (iterating over samples)...")
  for (j in seq_len(ncol(detP))) {
    message("  Processing sample ", j, " of ", ncol(detP), "...")

    # Type I Red
    intensity <- Red[TypeI.Red$AddressA, j] + Red[TypeI.Red$AddressB, j]
    detP[TypeI.Red$Name, j] <- pnorm(
      q = intensity,
      mean = 2 * rMu[j],
      sd = 2 * rSd[j],
      lower.tail = FALSE
    )

    # Type I Green
    intensity <- Green[TypeI.Green$AddressA, j] +
                Green[TypeI.Green$AddressB, j]
    detP[TypeI.Green$Name, j] <- pnorm(
      q = intensity,
      mean = 2 * gMu[j],
      sd = 2 * gSd[j],
      lower.tail = FALSE
    )

    # Type II
    intensity <- Red[TypeII$AddressA, j] + Green[TypeII$AddressA, j]
    detP[TypeII$Name, j] <- pnorm(
      q = intensity,
      mean = rMu[j] + gMu[j],
      sd = rSd[j] + gSd[j],
      lower.tail = FALSE
    )
  }

  message("Detection p-values calculated successfully.")


  PassProbeDF <- data.frame("Fraction" = rowSums(detP < as.numeric(opt$cutoff))/ncol(detP))
  PassProbeDF$Fraction <- round(PassProbeDF$Fraction, 3)
  PassProbeDF$ProbeIndex  <-seq_len(nrow(PassProbeDF))

  failed.probes <- data.frame("FailedProbes" = rownames(PassProbeDF)[PassProbeDF$Fraction < opt$threshold] )

    
  g <- ggplot(PassProbeDF, aes(x = ProbeIndex, y = Fraction)) +
    geom_hex(bins = 100) +   # adjust bins as needed
    labs(
      title = "Fraction of samples in which each probe passes the detection P-value cutoff",
      subtitle = paste0("Cutoff: ", opt$cutoff, " (red dashed line indicates threshold)\n",
                        "Probes passing the cutoff in at least ", opt$threshold * 100, "% of samples pass the filter"),
      x = "Probes",
      y = paste0("Fraction of samples in which detection P-value < ", opt$cutoff)
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank()
    ) +
    geom_hline(yintercept = opt$threshold, color = "red", linetype = "dashed") +
    annotate("text", x = max(PassProbeDF$ProbeIndex-10, na.rm = TRUE),
            y = opt$threshold,
            label = opt$cutoff,
            hjust = -0.1, vjust = -0.5, color = "red")

  # Open PDF for plotting
  pdf(opt$pdf_output, width = 10, height = 7)
  message("Generating detection p-value plots...")
  print(g)

  dev.off()
  
  message("Filtering probes based on detection p-values...")
  
  # Create a data frame for flagged samples and save to CSV
  if (dim(failed.probes)[1] > 0) {
    write.csv(failed.probes, file = opt$flagged, row.names = FALSE)
    message(paste0("Flagged ", dim(failed.probes)[1], " probes with p-value < ",opt$cutoff," in ",opt$threshold*100,"% samples saved to ", opt$flagged))
  } else {
    message("No probes flagged with p-value < ", opt$cutoff, " in ", opt$threshold * 100, "% of samples.")
    # Ensure an empty file is still created if no samples are flagged
    write.csv(data.frame("FailedProbes" = character()), file = opt$flagged, row.names = FALSE)
  }
  # Close the PDF device
  message("Detection p-value plots generated successfully and saved to ", opt$pdf_output)

}, error = function(e) {
  message("An error occurred: ", e$message)
})

