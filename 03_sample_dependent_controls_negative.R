#!/usr/bin/env Rscript
# =============================================================================
# 03_sample_dependent_controls_negative.R
# Author: Vartika Bisht
# Created: 13-Aug-2025
# Description:
# This script analyzes the performance of negative control probes to assess
# the system's background noise. It helps define the baseline signal floor,
# which is crucial for determining the detection limits of DNA methylation
# probes. This analysis is performed for both green and red channels.
# Arguments:
#    --rgset                → Path to the RGChannelSet object file (.RData)
#    --pdf_output           → Output PDF file path for the control probe plots
#    --platform             → DNA methylation platform (e.g., EPIC, 450K)
#    --manifestkey          → Path to the manifest key file (default: data/manifest.annotation.key.csv)
#    --assembly             → Assembly version of the DNA methylation platform
# Usage:
#   Rscript 03_sample_dependent_controls_negative.R --rgset <FILE.RData> \
#       --platform <PLATFORM> --assembly <ASSEMBLY> --pdf_output <output.pdf>
# Notes:
#   - Negative control probes are randomly permuted sequences not designed to
#     hybridize to DNA. Their average signal defines the system background, which
#     is a comprehensive measure of noise from various sources, including cross-hybridization
#     and imaging system background.
#   - The number of negative control probes varies by platform (e.g., 600 for
#     450K arrays and 407 for EPIC/EPICv2).
#   - While these probes are not used to flag samples directly, their analysis
#     is essential for a complete quality control assessment and for establishing
#     accurate detection limits for methylation probes.
# =============================================================================

## Usage
# ------

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
required_args <- c("rgset", "platform", "manifestkey", "assembly" , "pdf_output")

for (arg in required_args) {
  if (is.null(opt[[arg]])) {
    stop(paste("Error: Required argument --", arg, " is missing. Use --help for more information.", sep = ""))
  }
}

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

## Performance Metrics 
# -------------------
# | Purpose Name       | Count | Green Channel | Red Channel | Expected Value  |  Platform   |
# |--------------------|-------|---------------|-------------|-----------------|-------------|
# | Negative Average   | 600   |  (+)          |  (+)        | Background      | HM450K      |
# | Negative Average   | 407   |  (+)          |  (+)        | Background      | EPIC/EPICv2 |


ctrl = "NEGATIVE"
controlprobes = controlprobes[controlprobes$Type == ctrl &controlprobes$Color != -99,]
num_negative_probes = length(unique(controlprobes$ExtendedType))
message("Number of negative control probes: ", num_negative_probes)


red <- getRed(RGset)
green <- getGreen(RGset)
ctrlAddress <- getControlAddress(RGset, controlType = ctrl)
ctrlData <-  rbind(
            cbind(channel = "Red", melt(red[ctrlAddress, ], varnames = c("address", "sample"))),
            cbind(channel = "Green", melt(green[ctrlAddress, ], varnames = c("address", "sample"))))
ctrlData$value <- log2(ctrlData$value)
rm(red, green, ctrlAddress)

ctrlData <- ctrlData %>%
  group_by(sample, channel) %>%
  summarise(
    average = mean(value, na.rm = TRUE),
    std = sd(value, na.rm = TRUE),
    .groups = 'drop'
  )

ctrlData <- data.frame(ctrlData)
avg <- mean(ctrlData$average)

pdf(opt$pdf_output, width = 8, height = 11) 

ggplot(ctrlData, aes(x = channel, y = average, fill = channel)) +
  geom_boxplot(outlier.size = 0.5, alpha = 0.7) +
  scale_fill_manual(values = c("Red" = "red", "Green" = "green")) +
  labs(title = paste0("Average intensity of ",num_negative_probes," negative control probes per channel"),
       x = "",
       y = "Intensity") +
  theme_minimal() +
  scale_y_continuous(limits = c(6,15)) +
  theme(legend.position = "none",
  plot.title = element_text(hjust = 0.5)) + 
  geom_hline(yintercept = avg, linetype = "dashed", color = "black")

dev.off()
