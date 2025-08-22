#!/usr/bin/env Rscript
# =============================================================================
# 05_high_intensity_flagged_probes.R
# Author: Vartika Bisht & Tatiana Karp
# Created: 13-Aug-2025
# Description:
# This script identifies and flags probes with high signal intensity. Such probes
# are prone to signal saturation, which can lead to unreliable DNA methylation
# measurements, often resulting in beta values near 0.5. The script filters
# probes whose median intensity exceeds a specified cutoff.
# Arguments:
#    --mset                 → Path to the MethylSet object file (.RData)
#    --flagged              → Path to a CSV file to save the list of flagged probes
#    --cutoff               → Numeric threshold for filtering probes with high intensity (default: 10000)
#    --pdf_output           → Output PDF file path for plots
# Usage:
#   Rscript 05_high_intensity_flagged_probes.R --mset <MSet.RData> \
#     --flagged <high_intensity_flagged_probes.csv> --cutoff 10000
# Notes:
#   - High intensity probes are a source of technical artifact. When signal
#     intensities are saturated, the relationship between methylated and unmethylated
#     signals becomes unreliable.
#   - This issue is more prevalent with Infinium Type I probes.
#   - The script flags probes where the median of the total signal (methylated +
#     unmethylated) across all samples is greater than the specified cutoff.
#   - Removing these probes is a recommended step to improve the reliability of
#     downstream methylation analysis.
# =============================================================================

suppressPackageStartupMessages({
  tryCatch({
    library(minfi)
    library(optparse)
    library(ggplot2)
  }, error = function(e) {
    stop("Package loading failed: ", e$message, call. = FALSE)
  })
})

# Parse command line arguments
option_list <- list(
  make_option(c("-m", "--mset"), type = "character", default = NULL,
              help = "Path to the RGset object file"),
  make_option(c("-o", "--flagged"), type = "character", default = NULL,
              help = "Path to save the list of failed probes"),
  make_option(c("-c", "--cutoff"), type = "numeric", default = 10000,
              help = "High intensity cutoff threshold for filtering probes (default: 10000)"),
  make_option(c("-p", "--pdf_output"), type="character", default=NULL,
                help="Output PDF file path", metavar="FILE")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)


# VALIDATE REQUIRED ARGUMENTS
required_args <- c("mset", "cutoff", "flagged" , "pdf_output")

for (arg in required_args) {
  if (is.null(opt[[arg]])) {
    stop(paste("Error: Required argument --", arg, " is missing. Use --help for more information.", sep = ""))
  }
}


# Load the RGset object. This file is assumed to also contain the MSet object.
load(opt$mset)

## Extract probe information and intensities
# Get the manifest to identify Type I probes
manifest <- getManifest(MSet)
ProbeI <- getProbeInfo(manifest, type = "I")$Name

# Extract methylated and unmethylated signal intensities from the MethylSet
Meth <- getMeth(MSet)
Unmeth <- getUnmeth(MSet)

# Calculate the median methylated and unmethylated intensities for each probe
Methmed <- apply(Meth, 1, median)
Unmethmed <- apply(Unmeth, 1, median)

# Calculate the combined median intensity (geometric mean) for each probe
MUmed <- sqrt(Methmed * Unmethmed)

# Identify all probes with combined median intensity above the specified cutoff
hi.MU <- names(which(MUmed > opt$cutoff))

# ---
## Filter and save flagged probes

# Select only Type I probes from the high-intensity list
flagged_probes <- hi.MU[hi.MU %in% ProbeI]

# Save the list of flagged probes to the specified output file
write.csv(data.frame("FlaggedProbes" = flagged_probes), file = opt$flagged, quote = FALSE, row.names = FALSE)

print(paste("Number of probes flagged due to high intensity: ", length(flagged_probes)))
print(paste("List of flagged probes saved to: ", opt$flagged))


# Open PDF for plotting
pdf(opt$pdf_output, width = 10, height = 7)

# geometric mean intensities for each probe
MUmed <- data.frame(Probe = names(MUmed), GM = MUmed)

message("Generating intesity plot for type I probes...")
ggplot(MUmed, aes(x = GM)) +
  geom_histogram(bins = 50, fill = "blue", color = "black") +
  labs(title = "Distribution of Geometric Mean Intensity for Type I Probes",
        subtitle = paste0("Cutoff: ", opt$cutoff, " (red line indicates threshold)\n",
                          "Probes with geometric mean intensity above this cutoff are flagged"),
        x = "Geometric Mean Intensity",
        y = "Number of Probes") +
  theme_minimal() + 
  geom_vline(xintercept = opt$cutoff, color = "red", linetype = "dashed") +
  annotate("text", x = opt$cutoff, y = max(table(MUmed$GM)), 
            label = paste0("Cutoff ",opt$cutoff), color = "red", vjust = -1)

dev.off()

message("Filtering probes based on intensity...")
