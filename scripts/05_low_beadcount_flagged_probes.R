#!/usr/bin/env Rscript
# =============================================================================
# 05_low_beadcount_flagged_probes.R
# Author: Vartika Bisht & Tatiana Karp
# Created: 13-Aug-2025
# Description:
# This script identifies and flags probes with consistently low bead counts across
# samples. As Infinium array technology relies on multiple bead replicates per
# probe for accurate signal measurement, a low bead count (typically less than 3)
# indicates an unreliable signal that should be filtered out to improve data quality.
# Arguments:
#    --rgsetext             → Path to the ExtendedRGChannelSet object file (.RData)
#    --flagged              → Output CSV file to list flagged probes
#    --cutoff               → Percentage threshold for filtering probes (default: 5)
#    --pdf_output           → Output PDF file path for plots
# Usage:
#   Rscript 05_low_beadcount_flagged_probes.R --rgsetext <RGsetEXT.RData> \
#     --flagged <low_beadcount_flagged_probes.csv> --cutoff 5 \
#     --pdf_output <beadcount_plot.pdf>
# Notes:
#   - Illumina's Infinium arrays use bead replicates to measure the signal for each
#     CpG site. A sufficient number of bead replicates is crucial for obtaining
#     a reliable and precise measurement.
#   - This script flags probes where more than a specified percentage (cutoff)
#     of samples have a bead count less than 3.
#   - Filtering probes with low bead counts is an important quality control step
#     to remove unreliable data points that can affect downstream analysis.
# =============================================================================

suppressPackageStartupMessages({
  tryCatch({
    library(wateRmelon)
    library(minfi)
    library(optparse)
  }, error = function(e) {
    stop("Package loading failed: ", e$message, call. = FALSE)
  })
})


# Parse command line arguments

option_list <- list(
  make_option(c("-r", "--rgsetext"), type = "character", default = NULL,
              help = "Path to the RGset extended object file"),
  make_option(c("-o", "--flagged"), type = "character", default = NULL,
              help = "Path to save the list of failed probes"),
  make_option(c("-c", "--cutoff"), type = "numeric", default = 5,
              help = "Percentage threshold for filtering probes (default: 5)"),
  make_option(c("-p", "--pdf_output"), type = "character", default = NULL,
              help = "Path to save the histogram plot of average bead counts (PDF)")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)


# VALIDATE REQUIRED ARGUMENTS
required_args <- c("rgsetext", "cutoff", "flagged" , "pdf_output")

for (arg in required_args) {
  if (is.null(opt[[arg]])) {
    stop(paste("Error: Required argument --", arg, " is missing. Use --help for more information.", sep = ""))
  }
}

# Load the RGsetEXT object
load(opt$rgsetext)

## Perform bead count filtering
# Generate a matrix of bead counts. NAs represent probes with bead count < 3.
df.bead.counts <- wateRmelon::beadcount(RGsetEXT)

# Define a function to calculate the percentage of NAs in each row (probe)
pct.of.nas <- function(row) {
  sum(is.na(row)) / length(row) * 100
}

# Calculate the percentage of NAs for each probe
df.bead.counts$na_pct <- apply(df.bead.counts, 1, pct.of.nas)

# Identify probes that exceed the specified percentage of low bead counts
flagged_probes <- rownames(df.bead.counts[which(df.bead.counts$na_pct > opt$cutoff),])

# Save the list of flagged probes to the specified output file
write.csv(flagged_probes, file = opt$flagged, quote = FALSE, row.names = FALSE)

message(paste("Number of probes flagged due to low bead counts: ", length(flagged_probes)))
message(paste("List of flagged probes saved to: ", opt$flagged))

## Generate and save the histogram plot
# Calculate the average bead count for each probe, ignoring NA values.
# The `na.rm = TRUE` argument ensures that NAs (bead count < 3) are not included in the average.
if(length(flagged_probes) == 0) {
  message("No probes were flagged.")
  pdf(file = opt$pdf_output, width = 8, height = 6)
  dev.off()
} else{
  message("Generating histogram of average bead counts...")
  average_bead_counts <- rowMeans(df.bead.counts[, -ncol(df.bead.counts)], na.rm = TRUE)

  # Open a PDF device for plotting
  pdf(file = opt$pdf_output, width = 8, height = 6) # Set dimensions for the PDF

  # Create the histogram
  hist(average_bead_counts,
      main = "Histogram of Average Bead Counts Per Probe", # Main title of the plot
      xlab = "Average Number of Beads",                     # X-axis label
      ylab = "Occurrence (Number of Probes)",              # Y-axis label
      col = "skyblue",                                     # Bar color
      border = "black",                                    # Bar border color
      breaks = 50)                                         # Number of bins for the histogram

  # Close the PDF device, saving the plot
  dev.off()

  message(paste("Histogram plot saved to: ", opt$pdf_output))

}
