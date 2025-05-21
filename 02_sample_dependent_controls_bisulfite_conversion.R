#!/usr/bin/env Rscript
###############################################################################
# Rscript 02_sample_dependent_controls_bisulfite_conversion.R -r RGset.RData -c controls_probes.csv --bs1_output Flagged_bisulfite_conversion_I.csv --bs2_output Flagged_bisulfite_conversion_II.csv -p sample_dependent_controls_bisulfite_conversion_overview.pdf
#
# This script evaluates bisulfite conversion efficiency using Illumina's 
# sample-dependent control probes. The analysis is performed at both the
# sample and plate level to identify potential conversion issues.
#
# Methodology:
# 1. Uses control probes specifically designed to assess bisulfite conversion
#    efficiency (Type I and Type II)
# 2. Analyzes conversion performance per plate (since conversion is performed
#    at the plate level in experimental workflow)
# 3. Flags samples showing abnormal conversion patterns for further investigation
#
# Inputs:
# -r /path/to/RGset.RData          : Raw methylation data (RGChannelSet)
# -c /path/to/controls_probes.csv  : Control probe data
#
# Outputs:
# --bs1_output    : CSV of flagged samples from Type I probes
# --bs2_output    : CSV of flagged samples from Type II probes  
# -p              : PDF summary of conversion quality across plates
#
# References:
# Illumina Technical Note: 
# https://support-docs.illumina.com/ARR/Inf_HD_Methylation/Content/ARR/Methylation/SampleDependentControlsIntro_fINF_mMeth.htm
#
# Note: Bisulfite conversion occurs at the plate level prior to array hybridization.
###############################################################################


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
  make_option(c("-o", "--bs1_output"), type="character", default=NULL,
              help="Output file for BISULFITE CONVERSION I failures"),
  make_option(c("-t", "--bs2_output"), type="character", default=NULL,
              help="Output file for BISULFITE CONVERSION II failures"),
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

##############################################################
# Bisulfite Conversion Control Probes - Quality Assessment
###############################################################

# Purpose:
# Evaluate bisulfite conversion efficiency using Illumina's 
# control probe system. Two distinct probe types (I and II)
# assess conversion quality through different mechanisms.

# ==================== Key Concepts ====================

# Probe Design Principles:
# - Primer landing locations contain no cytosines (C)
# - Query regions contain interrogation sites
# - Type I: Highly sensitive to incomplete conversion
#   (Unconverted C hybridizes to M bead → false methylation signal)
# - Type II: Incomplete conversion affects single-base extension,
#   skewing beta values

# ==================== Control Probe System ====================

# Bisulfite Conversion I:
# - 12 control probes (6 methylated [C1-C6], 6 unmethylated [U1-U6])
# - Evaluates type I probe conversion efficiency

# Bisulfite Conversion II: 
# - 4 control probes
# - Evaluates type II probe conversion efficiency

# ==================== Type I Probe Specifications ====================

# Probe    Address    Color Channel   Type      Expected Signal
# ------------------------------------------------------------
# C1       22711390   Green           High      High (Methylated)
# C2       22795447   LimeGreen       High      High (Methylated) 
# C3       56682500   Lime            High      High (Methylated)
# U1       46651360   Blue            Background Low (Unmethylated)
# U2       24637490   SkyBlue         Background Low (Unmethylated)
# U3       33665449   Cyan            Background Low (Unmethylated)
# C4       54705438   Purple          High      High (Methylated)
# C5       49720470   Red             High      High (Methylated)
# C6       26725400   Tomato          High      High (Methylated)
# U4       57693375   Orange          Background Low (Unmethylated)
# U5       15700381   Gold            Background Low (Unmethylated)
# U6       33635504   Yellow          Background Low (Unmethylated)

# Channel Mapping:
# - C1-C3, U1-U3: Evaluated in Green channel
# - C4-C6, U4-U6: Evaluated in Red channel

# ==================== Quality Interpretation ====================

# Successful Conversion:
# - U1-U6: Low intensity (proper conversion of unmethylated Cs)
# - C1-C6: High intensity (retention of methylated Cs)

# Failed Conversion Indicators:
# - High intensity in U1-U6 probes OR
# - Low intensity in C1-C6 probes

# Technical Note:
# All color assignments in the manifest ultimately map to either
# green or red channels during signal processing.
##############################################################
ctrl = "BISULFITE CONVERSION I"
red <- getRed(RGset)
green <- getGreen(RGset)
ctrlAddress <- getControlAddress(RGset, controlType = ctrl)
ctrlData <-  rbind(
            cbind(channel = "Red", melt(red[ctrlAddress, ], varnames = c("address", "sample"))),
            cbind(channel = "Green", melt(green[ctrlAddress, ], varnames = c("address", "sample"))))
ctrlData$value <- log2(ctrlData$value)
ctrlData$sample_plate <- RGset@colData[ctrlData$sample,]$Sample_Plate
rm(red, green, ctrlAddress)
gc()

probe_groups <- control_probes[control_probes$V2 == ctrl,]
probe_groups$V4 <- substr(probe_groups$V4,17,19)
probe_groups$V5 <- ifelse(grepl("^C", probe_groups$V4), "Methylated", "Unmethylated")
colnames(probe_groups) <- c("address", "control_type", "probe_type", "color_label", "methylation_status")

# Merge with your intensity data
ctrlData <- ctrlData %>%
  left_join(probe_groups, by = "address")

# Check plate wise distribution of intensities
# Add high and background info with the expected channel of interest
ctrlData <- data.frame(ctrlData %>%
  mutate(  probe_group = case_when(
           color_label %in% c("C1", "C2", "C3") ~ "Green Channel High",
           color_label %in% c("U1", "U2", "U3") ~ "Green Channel Background",
           color_label %in% c("C4", "C5", "C6") ~ "Red Channel High",
           color_label %in% c("U4", "U5", "U6") ~ "Red Channel Background"
         )))

# Create the box plot to check the 
# Measure only in one channel, as specified before
df = ctrlData[ctrlData$channel == unlist(lapply(strsplit(ctrlData$probe_group, " "), function(x) x[1])),]
sample_dependent_controls_bisulfite_conversion_I_plate <- 
  ggplot(df, aes(x = sample_plate, y = value, fill = probe_group)) +
  geom_boxplot(position = position_dodge(0.8), width = 0.6) +
  facet_wrap(~sample_plate, scales = "free_x", ncol = 5) + 
  scale_fill_manual(values = c(
    "Green Channel High" = "darkgreen",
    "Green Channel Background" = "lightgreen",
    "Red Channel High" = "darkred",
    "Red Channel Background" = "pink"
  )) +
  labs(
    title = "Bisulfite Conversion Controls type I by plate",
    x = "Sample Plate",
    y = "Log2 Intensity",
    fill = "Probe Group"
  ) +
  theme_minimal() +
  theme(
    axis.title.x = element_blank(),     
    axis.text.x = element_blank(),      
    axis.ticks.x = element_blank(),   
    strip.background = element_blank(),  
    panel.spacing = unit(0.5, "lines")  
  )

rm(df)

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
sample_dependent_controls_bisulfite_conversion_I <- ggplot(ctrlData, aes(x = green, y = red, color = color_label)) +
  # Add diagonal x=y line first (so it's behind points)
  geom_abline(intercept = 0, slope = 1, color = "gray50", linetype = "dashed") +
  # Plot points
  geom_point(size = 1, alpha = 0.3) +
  # Make axes equal and set limits
  coord_equal(xlim = c(5,15), ylim = c(5, 15)) +
  labs(
    title = "Bisulfite Conversion Controls type I",
    subtitle = "Each sample is represented by 12 points (C1,C2,C3,C4,C5,C6,U1,U2,U3,U4,U5,U6)",
    x = "Log2 Green Channel Intensity",
    y = "Log2 Red Channel Intensity",
    color = "Probe Type"
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

# Flagging bad samples
# The C1-C6 should be high and U1-U6 should be low
# Calculate mean background intensity for controls
# U1-U3 represent background for green channel
bkg_green <- ctrlData %>%
  filter(color_label %in% c("U1", "U2", "U3")) %>%
  pull(green) %>%
  .[is.finite(.)] %>%
  mean(na.rm = TRUE)

# U4-U6 represent background for red channel
bkg_red <- ctrlData %>%
  filter(color_label %in% c("U4", "U5", "U6")) %>%
  pull(red) %>%
  .[is.finite(.)] %>%
  mean(na.rm = TRUE)

# Identify samples where controls (C1-C6) are below background thresholds
bad_samples <- ctrlData %>%
  group_by(sample) %>%
  summarize(
    green_low = mean(green[color_label %in% c("C1", "C2", "C3")], na.rm = TRUE) < bkg_green,
    red_low = mean(red[color_label %in% c("C4", "C5", "C6")], na.rm = TRUE) < bkg_red,
    .groups = "drop"
  ) %>%
  filter(green_low | red_low) %>%
  distinct(sample)

# Write output
write.csv(
  na.omit(bad_samples),
  file = opt$bs1_output,
  row.names = FALSE,
  col.names = FALSE
)


# ==================== Bisulfite Conversion II Controls ====================

# Purpose:
# Evaluate bisulfite conversion efficiency specifically for Type II probes
# using a set of 4 dedicated control probes.

# Key Characteristics:
# - All 4 probes are evaluated exclusively in the RED channel
# - Each probe has distinct color assignments in the manifest
# - Expected to show HIGH intensity signals when conversion is successful

# -------------------- Probe Specifications --------------------

# Probe ID    Address    Manifest Color    Probe Name
# --------------------------------------------------
# 1           43720395   Purple            BS Conversion II-1
# 2           70664314   Red               BS Conversion II-2
# 3           71718498   Tomato            BS Conversion II-3
# 4           30724412   Orange            BS Conversion II-4

# -------------------- Experimental Design --------------------

# Channel Evaluation:
# Probe Group    Green Channel    Red Channel    Expected Result
# -------------------------------------------------------------
# All 4 probes   - (Not used)     +              High Intensity

# Quality Interpretation:
# - Successful conversion: Consistently HIGH signals across all 4 probes
# - Potential issues: Low or variable signals among the probes

# Technical Note:
# Despite different color assignments in the manifest (Purple, Red, Tomato, Orange),
# all signals are ultimately measured through the RED channel during analysis.
ctrl = "BISULFITE CONVERSION II"
red <- getRed(RGset)
green <- getGreen(RGset)
ctrlAddress <- getControlAddress(RGset, controlType = ctrl)
ctrlData <-  rbind(
            cbind(channel = "Red", melt(red[ctrlAddress, ], varnames = c("address", "sample"))),
            cbind(channel = "Green", melt(green[ctrlAddress, ], varnames = c("address", "sample"))))
ctrlData$value <- log2(ctrlData$value)
ctrlData$sample_plate <- RGset@colData[ctrlData$sample,]$Sample_Plate

# Comparing intensity distribution platewise
# Everything is evaluated in red channel so remove green
df = ctrlData[ctrlData$channel == "Red",]

sample_dependent_controls_bisulfite_conversion_II_plate <- 
  ggplot(df, aes(y = sample_plate, x = value)) +  # Flip x and y
  geom_boxplot(width = 0.6,  fill = "lightblue") +  # Single color
  labs(
    title = "Bisulfite Conversion Controls type II by plate",
    y = "Sample Plate",  # Now on y-axis
    x = "Log2 Intensity red channel"  # Now on x-axis
  ) +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 8),  # Adjust plate label size
    panel.grid.major.y = element_blank(),  # Remove horizontal grid lines
    plot.title.position = "plot"  # Title alignment
  )

rm(df)



# Pivot to wide format (one row per address-sample pair)
ctrlData <- ctrlData %>%
  pivot_wider(
    names_from = channel,      # Column to split (Red/Green)
    values_from = value,       # Values to fill into new columns
    names_glue = "{tolower(channel)}"  # Force lowercase column names
  )

ctrlData <- data.frame(ctrlData)
rm(red, green, ctrlAddress)


# Get axis limits based on data range
axis_max <- max(c(ctrlData$green, ctrlData$red), na.rm = TRUE)
sample_dependent_controls_bisulfite_conversion_II <- ggplot(ctrlData, aes(x = green, y = red)) +
  # Add diagonal x=y line first (so it's behind points)
  geom_abline(intercept = 0, slope = 1, color = "gray50", linetype = "dashed") +
  # Plot points
  geom_point(size = 1, color = "gray80") +
  # Make axes equal and set limits
  coord_equal(xlim = c(5,15), ylim = c(5, 15)) +
  labs(
    title = "Bisulfite Conversion Controls type II",
    subtitle = "Each sample is represented by 4 points",
    x = "Log2 Green Channel Intensity",
    y = "Log2 Red Channel Intensity"
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

# Flag samples
# The 4 probes should be high
# Calculate mean intensity for all controls
# define background as mean - 2*sd
bkg = mean(ctrlData$red[is.finite(ctrlData$red)]) - 2*sd(ctrlData$red[is.finite(ctrlData$red)])

bad_samples = unique(ctrlData$sample)[unlist(lapply(unique(ctrlData$sample), function(x) {
  ifelse(mean(ctrlData[ctrlData$sample == x,]$red[is.finite(ctrlData[ctrlData$sample == x,]$red)]) < bkg, TRUE,FALSE)
}))]

# Write output
write.csv(
  na.omit(bad_samples),
  file = opt$bs2_output,
  row.names = FALSE,
  col.names = FALSE
)

# Finally we will check if there is any plate specifc bias in this step. It coudl be that all the samples we have flagged are actually coming from the same plate
pdf(opt$pdf_output, width = 8, height = 11)  # Increased height
grid.draw(ggplotGrob(sample_dependent_controls_bisulfite_conversion_I))
grid.newpage()
grid.draw(ggplotGrob(sample_dependent_controls_bisulfite_conversion_I_plate))
grid.newpage()
grid.draw(ggplotGrob(sample_dependent_controls_bisulfite_conversion_II))
grid.newpage()
grid.draw(ggplotGrob(sample_dependent_controls_bisulfite_conversion_II_plate))
dev.off()