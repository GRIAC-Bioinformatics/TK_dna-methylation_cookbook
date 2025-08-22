#!/usr/bin/env Rscript
# =============================================================================
# 03_sample_detection_p_value.R
# Author: Vartika Bisht
# Created: 13-Aug-2025
# Description:
# This script performs quality control by assessing detection p-values for each
# probe in a sample. It identifies samples with a high proportion of unreliable
# signals by comparing the total probe intensity to the background noise.
# Samples where a significant number of probes fail to meet a specified p-value
# cutoff are flagged for removal.
# Arguments:
#    --rgset                → Path to the RGChannelSet object file (.RData)
#    --detectionP           → Path to a file storing pre-calculated detection p-values (.RData)
#    --flagged              → Output CSV file to list flagged probes with high detection p-values
#    --cutoff               → Numeric cutoff for detection p-value (default: 0.01)
#    --threshold            → Threshold for the fraction of probes passing the p-value cutoff (default: 0.95)
#    --pdf_output           → Output PDF file path for QC plots
# Usage:
#   Rscript 03_sample_detection_p_value.R --rgset <FILE.RData> \
#     --detectionP <p_values.RData> --flagged <flagged_probes.csv> \
#     --cutoff 0.01 --threshold 0.95 --pdf_output <output.pdf>
# Notes:
#   - A detection p-value quantifies the reliability of a probe's signal. A low p-value
#     indicates the signal is significantly above background noise, which is calculated
#     using negative control probes.
#   - This script filters out samples where less than a specified fraction of probes
#     pass the detection p-value cutoff. For a cutoff of 0.01 and a threshold of 0.95,
#     a sample is flagged if more than 5% of its probes have a p-value greater than 0.01.
#   - The choice of cutoff and threshold values is empirical and should be tailored
#     to the specific dataset to balance data quality with sample retention.
# =============================================================================

# Silent library loading and error handling
suppressPackageStartupMessages({
  tryCatch({
  library(minfi)
  library(reshape2)
  library(optparse)
  library(ggplot2)
  }, error = function(e) {
    stop("Package loading failed: ", e$message, call. = FALSE)
  })
})

# Set up command line options
option_list <- list(
    make_option(c("-r", "--rgset"), type="character", default=NULL,
                help="Path to RGset file", metavar="FILE"),
    make_option(c("-d", "--detectionP"), type="character", default=NULL,
                help="Path to RData file storing Detection P values", metavar="FILE"),
    make_option(c("-f", "--flagged"), type="character", default=NULL,
                help="CSV file to save flagged probes which have high detection p-values", metavar="FILE"),
    make_option(c("-c", "--cutoff"), type="numeric", default=0.01,
                help="Cutoff for detection p-value (default: 0.01)", metavar="NUM"),
    make_option(c("-t", "--threshold"), type="numeric", default=0.95,
                help="Thershold for the fraction of probes which can fail per sample before the sample is flagged (default: 0.95)", metavar="NUM"),
    make_option(c("-p", "--pdf_output"), type="character", default=NULL,
                help="Output PDF file path", metavar="FILE")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)


# VALIDATE REQUIRED ARGUMENTS
required_args <- c("rgset", "detectionP", "flagged", "cutoff", "threshold" ,"pdf_output")

for (arg in required_args) {
  if (is.null(opt[[arg]])) {
    stop(paste("Error: Required argument --", arg, " is missing. Use --help for more information.", sep = ""))
  }
}


tryCatch({
  # Load the RGset object
  message("Loading RGset object...")
  load(opt$rgset) # Ensure opt$rgset is defined, e.g., from optparse
  message("=== Starting Detection P-value Calculation with minfi (custom) ===")

  message("Step 1: Extracting locus names...")
  locusNames <- getManifestInfo(RGset, "locusNames")
  message("   ✔ Extracted ", length(locusNames), " loci.")

  message("Step 2: Extracting negative control probe indices...")
  controlIdx <- getControlAddress(RGset, controlType = "NEGATIVE")
  message("   ✔ Found ", length(controlIdx), " negative controls.")

  message("Step 3: Extracting raw intensity signals (Red and Green channels)...")
  Red <- getRed(RGset)
  Green <- getGreen(RGset)
  message("   ✔ Red channel matrix: ", nrow(Red), " probes × ", ncol(Red), " samples.")
  message("   ✔ Green channel matrix: ", nrow(Green), " probes × ", ncol(Green), " samples.")

  message("Step 4: Extracting probe type information...")
  TypeI.Red <- getProbeInfo(RGset, type = "I-Red")
  TypeI.Green <- getProbeInfo(RGset, type = "I-Green")
  TypeII <- getProbeInfo(RGset, type = "II")
  message("   ✔ Type I-Red probes: ", nrow(TypeI.Red))
  message("   ✔ Type I-Green probes: ", nrow(TypeI.Green))
  message("   ✔ Type II probes: ", nrow(TypeII))

  message("Step 5: Setting up empty detection P-value matrix...")
  detP <- matrix(
      NA_real_,
      nrow = length(locusNames),
      ncol = ncol(Red),
      dimnames = list(locusNames, colnames(Red))
  )
  message("   ✔ Detection P-value matrix size: ", nrow(detP), " × ", ncol(detP))

  message("Step 6: Calculating background distributions (Red/Green from NEGATIVE controls)...")
  rBg <- Red[controlIdx, , drop = FALSE]
  gBg <- Green[controlIdx, , drop = FALSE]
  message("   ✔ Background extracted.")

  message("Step 7: Computing background medians and MADs...")
  rMu <- colMedians(rBg, na.rm = TRUE, useNames = TRUE)
  rSd <- matrixStats::colMads(rBg)
  gMu <- colMedians(gBg, na.rm = TRUE, useNames = TRUE)
  gSd <- matrixStats::colMads(gBg)
  message("   ✔ Red mean ± SD: ", signif(mean(rMu), 3), " ± ", signif(mean(rSd), 3))
  message("   ✔ Green mean ± SD: ", signif(mean(gMu), 3), " ± ", signif(mean(gSd), 3))

  message("Step 8: Filling detection P-value matrix sample by sample...")
  for (j in seq_len(ncol(detP))) {
      message("   → Processing sample ", j, "/", ncol(detP), " (", colnames(detP)[j], ") ...")

      # Type I Red
      intensity <- Red[TypeI.Red$AddressA, j] + Red[TypeI.Red$AddressB, j]
      detP[TypeI.Red$Name, j] <- pnorm(
          q = intensity,
          mean = 2 * rMu[j],
          sd = 2 * rSd[j],
          lower.tail = FALSE
      )

      # Type I Green
      intensity <- Green[TypeI.Green$AddressA, j] + Green[TypeI.Green$AddressB, j]
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

      message("      ✔ Completed sample ", j, ": ",
              sum(!is.na(detP[, j])), " probes assigned p-values.")
  }

  message("=== Detection P-value Calculation Completed Successfully ===")

  # Save the detection p-values matrix
  save(detP, file = opt$detectionP)

  fraction_of_probes_above_cutoff <- 1 - (colSums(detP>opt$cutoff)/nrow(detP))

  # Flag samples with mean detection p-value > opt$cutoff
  flagged_sample_names <- names(which(fraction_of_probes_above_cutoff < opt$threshold))
  num_flagged_samples <- length(flagged_sample_names)

  # Create a data frame for flagged samples and save to CSV
  if (num_flagged_samples > 0) {
    flagged_sample_df <- data.frame(Sample = flagged_sample_names)
    write.csv(flagged_sample_df, file = opt$flagged, row.names = FALSE)
    message(paste0("Flagged ", num_flagged_samples, " samples with less than ",100*opt$threshold, "% of it's probes passed with detection pvalue < ", opt$cutoff," saved to ", opt$flagged))
  } else {
    message("No samples flagged with mean p-value > ", opt$cutoff)
    # Ensure an empty file is still created if no samples are flagged
    write.csv(data.frame(Sample = character()), file = opt$flagged, row.names = FALSE)
  }


  # Open PDF for plotting
  pdf(opt$pdf_output, width = 10, height = 7)
    
  # 1. First, melt the Detection P value dataframe and apply the log10 transformation to the value column
  melted_detP <- melt(detP)
  melted_detP$value_log <- log10(melted_detP$value)

  # 2. Then, replace any -Inf values (which result from log10(0)) with -1000
  # The is.infinite() function is used to check for infinite values
  melted_detP$value_log[is.infinite(melted_detP$value_log)] <- -1000

  # Create the density plot using ggplot2
  # - The 'group' aesthetic is added to ensure a separate density curve is calculated for each sample.
  # - The 'color' is set to a static value outside of aes() to make all lines green.
  # - 'show.legend = FALSE' is added to hide the legend for Var2.
  print(ggplot(melted_detP, aes(x = value_log, group = Var2)) + 
    geom_density(color = "seagreen", show.legend = FALSE) +
    labs(
      title = "Density Plot of Log-Transformed Values",
      subtitle = paste0("Dashed red line represents the cutoff ",opt$cutoff),
      x = "Value (log10 scale)",
      y = "Density",
      # Add a caption to explain the data transformation
      caption = "Note: Values of 0, which result in a log-transformed value of -Inf, were changed to -1000."
    ) +
    theme_minimal() +
    # Add a slight padding to the bottom margin to ensure the caption is visible
    theme(plot.caption = element_text(size = 9, hjust = 0)) + 
    geom_vline(xintercept = opt$cutoff, color = "red" , linetype = "dashed"))


  dev.off()

  message("Detection p-value plots generated successfully and saved to ", opt$pdf_output)

}, error = function(e) {s
  message("An error occurred: ", e$message)
  if (!is.null(dev.list())) {
    dev.off()
  }
})