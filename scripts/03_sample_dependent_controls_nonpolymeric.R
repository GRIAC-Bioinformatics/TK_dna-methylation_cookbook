#!/usr/bin/env Rscript
# =============================================================================
# 03_sample_dependent_controls_nonpolymeric.R
# Author: Vartika Bisht
# Created: 13-Aug-2025
# Description:
# This script analyzes nonpolymorphic control probes to assess the overall
# performance of the DNA methylation assay from amplification to detection.
# It evaluates 4 control probes (one for each nucleotide: A, T, C, G) and
# identifies samples with low signal intensity, which may indicate technical issues.
# Arguments:
#    --rgset                → Path to the RGChannelSet object file (.RData)
#    --flagged              → Output CSV file path for flagged nonpolymorphic control failures
#    --pdf_output           → Output PDF file path for QC plots
#    --platform             → DNA methylation platform (e.g., EPIC, 450K)
#    --manifestkey          → Path to the manifest key file (default: data/manifest.annotation.key.csv)
#    --assembly             → Assembly version of the DNA methylation platform
# Usage:
#   Rscript 03_sample_dependent_controls_nonpolymeric.R \
#     --rgset <FILE.RData> --platform <PLATFORM> --assembly <ASSEMBLY> \
#     --flagged <flagged_nonpolymeric.csv> --pdf_output <output.pdf>
# Notes:
#   - Nonpolymorphic (NP) controls are designed to query a specific base (A, T, C, or G)
#     in a nonpolymorphic region of the bisulfite-converted genome.
#   - The four NP control probes are expected to have consistently high intensity signals.
#   - Channel-specific Evaluation:
#       - NP (A) and NP (T) probes are evaluated in the red channel.
#       - NP (C) and NP (G) probes are evaluated in the green channel.
#   - Quality Flagging: Samples are flagged as potentially low quality if the intensity
#     of any NP probe falls below a cutoff.
#   - Cutoff Definition: Since there is no recommended "low" signal value, the script
#     uses an empirical cutoff of mean intensity - (2 * standard deviation) for each probe.
# =============================================================================

# Silent library loading and error handling
suppressPackageStartupMessages({
  tryCatch({
    library(minfi)
    library(tidyverse)
    library(reshape2)
    library(dplyr) 
    library(ggplot2)
    library(optparse)
  }, error = function(e) {
    stop("Package loading failed: ", e$message, call. = FALSE)
  })
})


option_list <- list(
  make_option(c("-r", "--rgset"), type="character", default=NULL,
              help="Path to RGset RData file"),
  make_option(c("-f", "--flagged"), type="character", default=NULL,
              help="Output file for non polymeric controls failures"),
  make_option(c("-o", "--pdf_output"), type="character", default=NULL,
              help="Output PDF file for plots"),
  make_option(c("-p", "--platform"), type="character", default=NULL, 
              help="DNA methylation platform", metavar="character"),
  make_option(c("-k", "--manifestkey"), type="character", default="data/manifest.annotation.key.csv", 
              help="Manifest key file", metavar="character"),
  make_option(c("-a", "--assembly"), type="character", default=NULL,
              help="Assembly version of the DNA methylation platform", metavar="character")
)

# Parse command line options
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

# VALIDATE REQUIRED ARGUMENTS
required_args <- c("rgset", "platform", "manifestkey", "assembly" , "pdf_output", "flagged")

for (arg in required_args) {
  if (is.null(opt[[arg]])) {
    stop(paste("Error: Required argument --", arg, " is missing. Use --help for more information.", sep = ""))
  }
}

# Load the RGset data
message("Loading RGset data from: ", opt$rgset)
load(opt$rgset)



# LOAD MANIFEST AND ANNOTATION
message("Platform: ", opt$platform)
message("Assembly: ", opt$assembly)
message("Reading ",opt$manifestkey," ...")
manifest.key <- tryCatch({
  read.csv(opt$manifestkey, header = TRUE)
}, error = function(e) {
  stop("Failed to read manifest file: ", e$message)
})

Nmanifest <- manifest.key[manifest.key$platform == opt$platform & manifest.key$assembly == opt$assembly, ]$manifest
Nannotation <- manifest.key[manifest.key$platform == opt$platform & manifest.key$assembly == opt$assembly, ]$annotation

message("According to the ",opt$manifestkey," file, the following manifest and annotation will be used:")
message("Manifest: ", Nmanifest)
message("Annotation: ", Nannotation)

message("Loading the manifest and annotation packages...")
tryCatch({
    library(Nmanifest, character.only = TRUE)
    library(Nannotation, character.only = TRUE)
    data(list = Nannotation)
}, error = function(e) {
    stop("Failed to load manifest or annotation package: ", e$message)
})

message("Manifest and annotation packages loaded successfully.")
message("Getting control probe information from ",Nannotation)
controlprobes = getProbeInfo(get(Nannotation) , type = "Control")

# Load the RGset data
message("Loading RGset data from: ", opt$rgset)
load(opt$rgset)

# Nonpolymorphic (NP) controls: 450K
# Probe ID    Control Type    Color    Description
# 24701411    NON-POLYMORPHIC    Red    NP (A)
# 18773482    NON-POLYMORPHIC    Purple    NP (T)
# 23663352    NON-POLYMORPHIC    Green    NP (C)
# 70645401    NON-POLYMORPHIC    Blue    NP (G)

# Purpose Name on the Array    Number    Evaluate Green (GRN)    Evaluate Red (RED)    Expected Intensity
# Non-Polymorphic NP (A), (T)    2    -    +    High
# Non-Polymorphic NP (C), (G)    2    +    -    High

ctrl = "NON-POLYMORPHIC"
red <- getRed(RGset)
green <- getGreen(RGset)
ctrlAddress <- getControlAddress(RGset, controlType = ctrl)
ctrlData <-  rbind(
            cbind(channel = "Red", melt(red[ctrlAddress, ], varnames = c("address", "sample"))),
            cbind(channel = "Green", melt(green[ctrlAddress, ], varnames = c("address", "sample"))))
ctrlData$value <- log2(ctrlData$value)
rm(red, green, ctrlAddress)

probe_groups = controlprobes[controlprobes$Type == ctrl &controlprobes$Color != -99,]
probe_groups$ExtendedType <- substr(probe_groups$ExtendedType,5,5)
colnames(probe_groups) <- c("address", "control_type", "color_label", "base")
probe_groups$expected_channel = ifelse(probe_groups$base %in% c("A","T"),"Red","Green")
probe_groups <- data.frame(probe_groups)
probe_groups$address <- as.numeric(probe_groups$address)

# Merge with your intensity data
ctrlData <- ctrlData %>%
  left_join(probe_groups, by = "address")


# Only maintain the intensity values from the channels that are expected by Illumina
ctrlData <- ctrlData %>%
  filter(channel == expected_channel) 

ctrlData <- data.frame(ctrlData)


pdf(opt$pdf_output, width = 8, height = 11)
ggplot(ctrlData, aes(x = base, y = value, fill = expected_channel)) +
  geom_boxplot(
    outlier.color = "black",        # Make outliers clearly visible in red
    outlier.shape = 16,           # Solid circle for outliers
    outlier.size = 1.5,           # Slightly larger outlier points
    alpha = 0.8                   # Slight transparency for boxplots
  ) +
  scale_fill_manual(
    values = c("Red" = "firebrick2", "Green" = "seagreen3"),  # Color mapping
    name = "Expected Channel"      # Legend title
  ) +
  labs(
    title = "Non Polymeric Controls",
    x = "",
    y = "Log2 Intensity",
    subtitle = "Colored by expected channel"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5),
    legend.position = "bottom",
    panel.grid.major.x = element_blank()  # Remove vertical grid lines for cleaner look
  ) +
  scale_x_discrete(limits = c("A", "C", "T", "G"))  # Ensure all bases are shown in order

dev.off()


# Flag samples where signal is lower than background threshold 
# Calculate background threshold (mean - 2*SD of finite values)
bkg = mean(ctrlData$value[is.finite(ctrlData$value)]) - 2*sd(ctrlData$value[is.finite(ctrlData$value)])
# Identify bad samples below threshold
bad_samples = unique(ctrlData[ctrlData$value < bkg,]$sample)

# Write output
write.csv(
  bad_samples,
  file = opt$flagged,
  row.names = FALSE
)
