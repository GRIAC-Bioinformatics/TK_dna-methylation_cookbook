#Detection p-value filtering
# The same as for samples, we want to remove probes showing high detection p-values. 
# For example the scanner can failed to read signals from some probes due to low intensities or artifacts in the array.
# Filter probes with detection p-values higher than the threshold (cut off could be 0.01, but that is also arbitrary):
# It could be probes that have failed in one or more samples/ or probes that have failed in at least x% of the samples
# In this case we remove probes thathave detection p-value>0.01 in more than 1% of the samples: 
# Usage : Rscript 04_detection_pvalue_flagged_probes.R --rgset RGset.RData --flagged Flagged_probe_detectionP.csv --cutoff 0.01 --threshold 0.99 --pdf detectionP_flagged_probes.pdf


# Silent library loading and error handling
suppressPackageStartupMessages({
  tryCatch({
  library(minfi)
  library(optparse)
  library(ggplot2)
  library(matrixStats)
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


load(opt$rgset)
# This RGset object has been filtered for "bad samples" in the previous step
# RGset --> this is expected to be the name of the loaded object
message("Loading sample filtered RGset object...")

tryCatch({

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

  PassProbeDF <- data.frame("Fraction" = rowSums(detP < opt$cutoff)/ncol(detP))
  failed.probes <- data.frame("FailedProbes" = names(which(rowSums(detP < opt$cutoff)/ncol(detP) < opt$threshold)))

  # Open PDF for plotting
  pdf(opt$pdf, width = 10, height = 7)
  message("Generating detection p-value plots...")
  print(ggplot(PassProbeDF, aes(x = Fraction)) +
      geom_histogram(binwidth = 0.001, fill = "blue", color = "black") +
      labs(title = "Distribution of number of samples in which each probes passes the detection P-value cutoff",
          subtitle = paste0("Cutoff: ", opt$cutoff, " (red line indicates threshold)\n",
                            "Probes passing the cutoff in at least ", opt$threshold * 100, "% of samples are retained"),
          x = "Fraction of Samples Passing Cutoff",
          y = "Number of Probes") +
      theme_minimal() +
      scale_x_log10() + # Added log scale for x-axis
      geom_vline(xintercept = opt$threshold, color = "red", linetype = "dashed") +
      annotate("text", x = opt$threshold, y = max(table(PassProbeDF$Fraction)),
              label = paste0("Cutoff ",opt$threshold), color = "red", vjust = -1))

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
  message("Detection p-value plots generated successfully and saved to ", opt$pdf)

}, error = function(e) {
  message("An error occurred: ", e$message)
  # Ensure dev.off() is called if pdf() was opened to avoid corrupted files
  if (!is.null(dev.list())) {
    dev.off()
  }
})

