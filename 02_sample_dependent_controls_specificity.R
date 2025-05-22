#!/usr/bin/env Rscript
# Author: Vartika Bisht
# Date: 22 May 2025
# Description: 
#   This script analyzes specificity controls from methylation array data to monitor
#   non-specific primer extension for probe types I and II. Specificity controls
#   evaluate proper base incorporation during the assay process.

# Control Probe Details:
# ======================
# * Specificity I Controls (12 probes - measured in both Green and Red channels):
#   - PM (Perfect Match) probes end with A (A/T match - should have HIGH signal)
#   - MM (Mismatch) probes end with G (G/T mismatch - should have LOW signal)
#   - PM1-3/MM1-3 measured in Red channel
#   - PM4-6/MM4-6 measured in Green channel
#
# * Specificity II Controls (3 probes - measured in Red channel only):
#   - Should incorporate "A" base across nonpolymorphic T (Red channel signal)
#   - Non-specific "G" base incorporation reduces Red channel signal

# Purpose:
# ========
# These controls monitor primer annealing specificity after samples are placed on slides
# (not plate-specific). The analysis flags samples where:
# 1. PM signals are lower than MM signals (indicating specificity issues)
# 2. Red channel signals for Type II probes are abnormally low

# Usage:
# ======
# Rscript 02_sample_dependent_controls_specificity.R 
#   -r RGset.RData 
#   -c controls_probes.csv 
#   --sp1_output Flagged_specificity_I.csv 
#   --sp2_output Flagged_specificity_II.csv 
#   -p sample_dependent_controls_specificity_overview.pdf

# Parameters:
# ===========
# Input:
#   -r, --rgset     : Path to RGset.RData (Raw intensity data)
#   -c, --controls  : Path to controls_probes.csv (Control probe metadata)
#
# Output:
#   --sp1_output    : CSV file for flagged Type I specificity issues
#   --sp2_output    : CSV file for flagged Type II specificity issues
#   -p, --plot      : PDF report with control visualizations

# Note: 
# =====
# This analysis occurs after sample placement on slides and evaluates primer
# annealing specificity, not plate effects.

# Silent library loading and error handling
suppressPackageStartupMessages({
  tryCatch({
    library(minfi)
    library(IlluminaHumanMethylation450kmanifest)
    library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
    library(reshape2)
    library(ggplot2)
    library(dplyr)
    library(tidyr)
    library(grid)
    library(gridExtra)
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
  make_option(c("-o", "--sp1_output"), type="character", default=NULL,
              help="Output file for SPECIFICITY I failures"),
  make_option(c("-t", "--sp2_output"), type="character", default=NULL,
              help="Output file for SPECIFICITY II failures"),
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

# ==================== SPECIFICITY I Controls ====================

ctrl = "SPECIFICITY I"
red <- getRed(RGset)
green <- getGreen(RGset)
ctrlAddress <- getControlAddress(RGset, controlType = ctrl)
ctrlData <-  rbind(
            cbind(channel = "Red", melt(red[ctrlAddress, ], varnames = c("address", "sample"))),
            cbind(channel = "Green", melt(green[ctrlAddress, ], varnames = c("address", "sample"))))
ctrlData$value <- log2(ctrlData$value)
rm(red, green, ctrlAddress)

# Perfect match (PM) and mismatch (MM) probes
# PM should be higher than MM
probe_groups <- control_probes[control_probes$V2 == ctrl,]
probe_groups$V5 <- unlist(lapply(strsplit(probe_groups$V4, " "), function(x) x[3]))
probe_groups$V4 <- substr(probe_groups$V4,16,17)
probe_groups$V6 <- paste(probe_groups$V4,probe_groups$V5, sep = "")
colnames(probe_groups) <- c("address", "control_type", "color_label", "match_status", "group", "match")


# Merge with your intensity data
ctrlData <- ctrlData %>%
  left_join(probe_groups, by = "address")

# Pivot to wide format (one row per address-sample pair)
ctrlData <- ctrlData %>%
  pivot_wider(
    names_from = channel,      # Column to split (Red/Green)
    values_from = value,       # Values to fill into new columns
    names_glue = "{tolower(channel)}"  # Force lowercase column names
  )

ctrlData <- data.frame(ctrlData)

# Get axis limits based on data range
axis_max <- max(c(ctrlData$green, ctrlData$red), na.rm = TRUE)

sample_dependent_controls_specificity_I <- ggplot(ctrlData, aes(x = green, y = red, color = match)) +
  # Add diagonal x=y line first (so it's behind points)
  geom_abline(intercept = 0, slope = 1, color = "gray50", linetype = "dashed") +
  # Plot points
  geom_point(size = 1, alpha = 0.3) +
  # Make axes equal and set limits
  coord_equal(xlim = c(5,15), ylim = c(5, 15)) +
  labs(
    title = "Specificity Controls type I",
    subtitle = "Controls to be evaluated in both channels\nEach sample represented by 12 points\n(PM1, PM2, PM3, MM1, MM2, MM3 in red and PM4, PM5, PM6, MM4, MM5, MM6 in green)",
    x = "Log2 Green Channel Intensity",
    y = "Log2 Red Channel Intensity",
    color = "Match Type"
  ) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    legend.title = element_text(size = 10),
    plot.title = element_text(hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    aspect.ratio = 1  # Ensure square plot area
  ) +
  guides(color = guide_legend(override.aes = list(size = 3, alpha = 1)))

# Flag samples where PM signal is lower than MM signal
bad_samples <- ctrlData %>%
  group_by(sample) %>%
  summarize(
    # Calculate mean signals for each probe type
    green_PM_mean = mean(green[match %in% c("PM1", "PM2", "PM3")], na.rm = TRUE),
    green_MM_mean = mean(green[match %in% c("MM1", "MM2", "MM3")], na.rm = TRUE),
    red_PM_mean = mean(red[match %in% c("PM4", "PM5", "PM6")], na.rm = TRUE),
    red_MM_mean = mean(red[match %in% c("MM4", "MM5", "MM6")], na.rm = TRUE),
    
    # Flag if either channel has PM < MM
    is_bad = any(green_PM_mean < green_MM_mean, red_PM_mean < red_MM_mean)
  ) %>%
  filter(is_bad) %>%
  pull(sample) %>%
  na.omit()

# Write output
write.csv(
  bad_samples,
  file = opt$sp1_output,
  row.names = FALSE,
  col.names = FALSE
)

# ==================== SPECIFICITY II Controls ====================

ctrl = "SPECIFICITY II"
red <- getRed(RGset)
green <- getGreen(RGset)
ctrlAddress <- getControlAddress(RGset, controlType = ctrl)
ctrlData <-  rbind(
            cbind(channel = "Red", melt(red[ctrlAddress, ], varnames = c("address", "sample"))),
            cbind(channel = "Green", melt(green[ctrlAddress, ], varnames = c("address", "sample"))))
ctrlData$value <- log2(ctrlData$value)
rm(red, green, ctrlAddress)


# Merge with your intensity data
ctrlData <- ctrlData %>%
  left_join(probe_groups, by = "address")

# Pivot to wide format (one row per address-sample pair)
ctrlData <- ctrlData %>%
  pivot_wider(
    names_from = channel,      # Column to split (Red/Green)
    values_from = value,       # Values to fill into new columns
    names_glue = "{tolower(channel)}"  # Force lowercase column names
  )

ctrlData <- data.frame(ctrlData)

# The signal is suppose to be recorded only in the red channel
sample_dependent_controls_specificity_II <- ggplot(ctrlData, aes(x = red)) +
  geom_histogram(
    binwidth = 0.1,
    fill = "skyblue",
    color = "white",
    alpha = 0.8
  ) +
  labs(
    title = "Specificity Controls type II",
    subtitle = "Controls to be evaluated in the red channel\nEach sample represented by 3 points",
    x = "Log2 Red Channel Intensity",
    y = "Count"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5),
    panel.grid.major.x = element_line(color = "gray90"),
    panel.grid.minor.x = element_blank(),
    plot.subtitle = element_text(hjust = 0.5),
  )

# Flag samples where signal is lower than background threshold 
# Calculate background threshold (mean - 2*SD of finite red channel values)
bkg = mean(ctrlData$red[is.finite(ctrlData$red)]) - 2*sd(ctrlData$red[is.finite(ctrlData$red)])
# Identify bad samples below threshold
bad_samples = unique(ctrlData[ctrlData$red < bkg,]$sample)

# Write output
write.csv(
  bad_samples,
  file = opt$sp2_output,
  row.names = FALSE,
  col.names = FALSE
)

# Finally we will check if there is any plate specifc bias in this step. It coudl be that all the samples we have flagged are actually coming from the same plate
pdf(opt$pdf_output, width = 8, height = 11)  # Increased height
grid.draw(ggplotGrob(sample_dependent_controls_specificity_I))
grid.newpage()
grid.draw(ggplotGrob(sample_dependent_controls_specificity_II))
dev.off()
