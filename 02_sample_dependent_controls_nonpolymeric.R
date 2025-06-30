#!/usr/bin/env Rscript
# Title: Sample Dependent Controls - Nonpolymorphic (NP) Controls Analysis
# Author: Vartika Bisht
# Date: 22 May 2025
# Description: 
#   This script analyzes nonpolymorphic control probes from methylation array data.
#   Nonpolymorphic controls test the overall performance of the assay by querying
#   specific bases in nonpolymorphic regions of the bisulfite genome.
#   The script evaluates 4 control probes (one for each nucleotide: A,T,C,G)
#   and identifies samples with potential technical issues.
#
# Input Parameters:
#   -r, --rgset        : Path to RGset RData file containing raw intensity data
#   -c, --control_probes : Path to CSV file with control probe information
#   -o, --np_output    : Output CSV file path for flagged nonpolymorphic control failures
#   -p, --pdf_output   : Output PDF file path for control probe visualization
#
# Output:
#   - CSV file containing samples that failed nonpolymorphic control checks
#   - PDF report with quality control plots showing control probe intensities
#
# Usage Example:
#   Rscript 02_sample_dependent_controls_nonpolymeric.R \
#     -r RGset.RData \
#     -c controls_probes.csv \
#     -o Flagged_nonpolymeric.csv \
#     -p sample_dependent_controls_nonpolymeric_overview.pdf

# Nonpolymorphic (NP) controls explanation:
# Nonpolymorphic controls test the overall performance of the assay, from amplification 
# to detection, by querying a particular base in a nonpolymorphic region of the bisulfite genome.
# One non-polymorphic control has been designed to query each of the 4 nucleotides (A,T, C, and G)
# 4 non-polymorphic control probes, one for each nucleotide.

# Probe ID    Control Type    Color    Description
# 24701411    NON-POLYMORPHIC    Red    NP (A)
# 18773482    NON-POLYMORPHIC    Purple    NP (T)
# 23663352    NON-POLYMORPHIC    Green    NP (C)
# 70645401    NON-POLYMORPHIC    Blue    NP (G)

# Purpose Name on the Array    Number    Evaluate Green (GRN)    Evaluate Red (RED)    Expected Intensity
# Non-Polymorphic NP (A), (T)    2    -    +    High
# Non-Polymorphic NP (C), (G)    2    +    -    High

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
  make_option(c("-c", "--control_probes"), type="character", default=NULL,
              help="Path to Illumina control probes CSV file"),
  make_option(c("-o", "--np_output"), type="character", default=NULL,
              help="Output file for non polymeric controls failures"),
  make_option(c("-p", "--pdf_output"), type="character", default=NULL,
              help="Output PDF file for plots")
)

# Parse command line options
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

# Load the RGset data
message("Loading RGset data from: ", opt$rgset)
load(opt$rgset)

# Check if RGset exists after loading
if (!exists("RGset")) {
  stop("RGset object not found in the RData file. Please check the input file.")
}

control_probes <- read.csv(opt$control_probes, header = FALSE)[,c(1,2,3,4)]


ctrl = "NON-POLYMORPHIC"
red <- getRed(RGset)
green <- getGreen(RGset)
ctrlAddress <- getControlAddress(RGset, controlType = ctrl)
ctrlData <-  rbind(
            cbind(channel = "Red", melt(red[ctrlAddress, ], varnames = c("address", "sample"))),
            cbind(channel = "Green", melt(green[ctrlAddress, ], varnames = c("address", "sample"))))
ctrlData$value <- log2(ctrlData$value)
rm(red, green, ctrlAddress)

probe_groups <- control_probes[control_probes$V2 == ctrl,]
probe_groups$V4 <- substr(probe_groups$V4,5,5)
colnames(probe_groups) <- c("address", "control_type", "color_label", "base")
probe_groups$expected_channel = ifelse(probe_groups$base %in% c("A","T"),"Red","Green")

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
  file = opt$np_output,
  row.names = FALSE,
  col.names = FALSE
)
