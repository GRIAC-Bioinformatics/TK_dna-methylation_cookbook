#!/usr/bin/env Rscript
# =============================================================================
# 03_sample_dependent_controls_bisulfite_conversion.R
# Author: Vartika Bisht
# Created: 13-Aug-2025
# Description:
# This script evaluates the efficiency of bisulfite conversion using Illumina's
# sample-dependent control probes. It generates plots and data tables to
# assess conversion quality at both the sample and plate levels, which is
# critical for identifying technical issues in DNA methylation experiments.
# Arguments:
#    --rgset                → Path to the RGChannelSet object file (.RData)
#    --platform             → DNA methylation platform (e.g., EPIC, 450K)
#    --manifestkey          → Path to the manifest key file (default: data/manifest.annotation.key.csv)
#    --assembly             → Assembly version of the DNA methylation platform
#    --bs1_output           → Output CSV file path for Bisulfite Conversion I failures
#    --bs2_output           → Output CSV file path for Bisulfite Conversion II failures
#    --pdf_output           → Output PDF file path for QC plots
# Usage:
#   Rscript 03_sample_dependent_controls_bisulfite_conversion.R \
#       --rgset <FILE.RData> --platform <PLATFORM> --assembly <ASSEMBLY> \
#       --bs1_output <bs1_failures.csv> --bs2_output <bs2_failures.csv> \
#       --pdf_output <output.pdf>
# Notes:
#   - This script assesses the efficiency of bisulfite conversion on genomic DNA
#     using a set of synthesized control probes that undergo the same conversion step.
#   - Bisulfite Conversion I Probes:
#      - These probes are designed to have high intensity in the C probes (C1-C5) and low intensity in the U probes (U1-U5).
#      - Channel Mapping: Theoretical colors in the manifest are mapped to two channels:
#        - Green Channel: C1, C2, U1, and U2.
#        - Red Channel: C3, C4, C5, U3, U4, and U5.
#      - Flagging Condition: A sample is flagged if its mean C probe intensity is less than its mean U probe intensity in either the green or red channel.
#   - Bisulfite Conversion II Probes:
#      - There are 4 Type II control probes used to check bisulfite conversion efficiency.
#      - All 4 probes are evaluated in the red channel and are expected to show high intensity signals.
#      - Since there is no absolute definition of a "high" signal, this script defines it as a range based on the mean and standard deviation of the probe intensities across all samples.
# =============================================================================

# Silent library loading and error handling
suppressPackageStartupMessages({
  tryCatch({
    library(minfi)
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
  make_option(c("-p", "--platform"), type="character", default=NULL, 
              help="DNA methylation platform", metavar="character"),
  make_option(c("-k", "--manifestkey"), type="character", default="data/manifest.annotation.key.csv", 
              help="Manifest key file", metavar="character"),
  make_option(c("-a", "--assembly"), type="character", default=NULL,
              help="Assembly version of the DNA methylation platform", metavar="character"),
  make_option(c("-o", "--bs1_output"), type="character", default=NULL,
              help="Output file for BISULFITE CONVERSION I failures"),
  make_option(c("-t", "--bs2_output"), type="character", default=NULL,
              help="Output file for BISULFITE CONVERSION II failures"),
  make_option(c("-f", "--pdf_output"), type="character", default=NULL,
              help="Output PDF file for plots")
)

# Parse command line options
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

# VALIDATE REQUIRED ARGUMENTS
required_args <- c("rgset", "platform", "manifestkey", "assembly" , "bs1_output" , "bs2_output", "pdf_output")

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

green_channel_label = c("Green","LimeGreen","Lime","Blue","SkyBlue","Cyan")
red_channel_label = c("Purple","Red","Tomato","Orange","Gold","Yellow")
message("Assuming the following channel labels:")
message("Theoretical colors measured in green channel : " ,paste0(green_channel_label, collapse = " , "))
message("Theoretical colors measured in red channel : " ,paste0(red_channel_label, collapse = " , "))

message("Loading RGChannelSet object from ", opt$rgset)
load(opt$rgset)

# =============================================================
# Bisulfite Conversion Control Probes - Quality Assessment
# =============================================================

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
# https://support-docs.illumina.com/ARR/Inf_HTS_Methylation/Content/ARR/Methylation/ControlBeadTypeIDs_fINF_mMeth.htm#BeadTypeIDs
# EPIC arrays, C1, C2 are green channel controls & C3-C5 are treated as red channel controls in "Bisulfite Conversion I" in the Controls table: https://help.dragenarray.illumina.com/dragen-array-v1.0/product-guides/output-files
# Illumina dropped U6 & C6 when they made the EPIC v1.0 array; see "Removed one pair of low intensity bisulfite controls" here: https://support.illumina.com/content/dam/illumina-support/documents/downloads/productfiles/methylationepic/infinium-methylationepic-v1-0-b5-customer-release-notes.pdf
# 450K (Infinium HumanMethylation450 BeadChip)
# Address    Type                     Color       ExtendedType
# -------------------------------------------------------------
# 22711390   BISULFITE CONVERSION I   Green       BS Conversion I-C1
# 22795447   BISULFITE CONVERSION I   LimeGreen   BS Conversion I-C2
# 56682500   BISULFITE CONVERSION I   Lime        BS Conversion I-C3
# 54705438   BISULFITE CONVERSION I   Purple      BS Conversion I-C4
# 49720470   BISULFITE CONVERSION I   Red         BS Conversion I-C5
# 26725400   BISULFITE CONVERSION I   Tomato      BS Conversion I-C6
# 46651360   BISULFITE CONVERSION I   Blue        BS Conversion I-U1
# 24637490   BISULFITE CONVERSION I   SkyBlue     BS Conversion I-U2
# 33665449   BISULFITE CONVERSION I   Cyan        BS Conversion I-U3
# 57693375   BISULFITE CONVERSION I   Orange      BS Conversion I-U4
# 15700381   BISULFITE CONVERSION I   Gold        BS Conversion I-U5
# 33635504   BISULFITE CONVERSION I   Yellow      BS Conversion I-U6


# EPIC (Infinium MethylationEPIC 850K BeadChip)
# Address    Type                     Color     ExtendedType
# ------------------------------------------------------------
# 22795447   BISULFITE CONVERSION I   Green     BS Conversion I-C1
# 56682500   BISULFITE CONVERSION I   Lime      BS Conversion I-C2
# 54705438   BISULFITE CONVERSION I   Purple    BS Conversion I-C3
# 49720470   BISULFITE CONVERSION I   Red       BS Conversion I-C4
# 26725400   BISULFITE CONVERSION I   Tomato    BS Conversion I-C5
# 24637490   BISULFITE CONVERSION I   Blue      BS Conversion I-U1
# 33665449   BISULFITE CONVERSION I   Cyan      BS Conversion I-U2
# 57693375   BISULFITE CONVERSION I   Orange    BS Conversion I-U3
# 15700381   BISULFITE CONVERSION I   Gold      BS Conversion I-U4
# 33635504   BISULFITE CONVERSION I   Yellow    BS Conversion I-U5


# EPICv2 (Infinium MethylationEPIC v2 950K BeadChip)
# Address    Type                     Color     ExtendedType
# ------------------------------------------------------------
# 22795447   BISULFITE CONVERSION I   Green     BS Conversion I-C1
# 56682500   BISULFITE CONVERSION I   Lime      BS Conversion I-C2
# 54705438   BISULFITE CONVERSION I   Purple    BS Conversion I-C3
# 49720470   BISULFITE CONVERSION I   Red       BS Conversion I-C4
# 26725400   BISULFITE CONVERSION I   Tomato    BS Conversion I-C5
# 24637490   BISULFITE CONVERSION I   Blue      BS Conversion I-U1
# 33665449   BISULFITE CONVERSION I   Cyan      BS Conversion I-U2
# 57693375   BISULFITE CONVERSION I   Orange    BS Conversion I-U3
# 15700381   BISULFITE CONVERSION I   Gold      BS Conversion I-U4
# 33635504   BISULFITE CONVERSION I   Yellow    BS Conversion I-U5
#
# Channel Mapping:
# - C1-C3, U1-U3: Evaluated in Green channel
# - C4-C6, U4-U6: Evaluated in Red channel
# However when it comes to EPIC and EPICv2 :
# Channel Mapping:
# - C1-C2, U1-U2: Evaluated in Green channel
# - C3-C5, U3-U5: Evaluated in Red channel
# Channel expectation based on internal color label
# Assuming :
# Green channel color : Green, LimeGreen , Lime , Blue , SkyBlue , Cyan
# Red channel color : Purple , Red , Tomato , Orange , Gold , Yellow

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

ctrl = "BISULFITE CONVERSION I"
red <- getRed(RGset)
green <- getGreen(RGset)
ctrlAddress <- getControlAddress(RGset, controlType = ctrl)
ctrlData <-  rbind(
            cbind(channel = "Red", melt(red[ctrlAddress, ], varnames = c("address", "sample"))),
            cbind(channel = "Green", melt(green[ctrlAddress, ], varnames = c("address", "sample"))))
ctrlData$value <- log2(ctrlData$value)
ctrlData$sample_plate <- RGset@colData[ctrlData$sample,]$Plate
ctrlData$sample_well <- RGset@colData[ctrlData$sample,]$Well
ctrlData$slide <- RGset@colData[ctrlData$sample,]$Slide
ctrlData$array <- RGset@colData[ctrlData$sample,]$Array

probe_groups <- controlprobes[controlprobes$Type == ctrl,]
probe_groups$ExtendedType <- substr(probe_groups$ExtendedType,17,19)
probe_groups$methylation_status <- ifelse(grepl("^C", probe_groups$ExtendedType), "Methylated", "Unmethylated")
colnames(probe_groups) <- c("address", "control_type", "color_label", "probe_type", "methylation_status")
probe_groups <- data.frame(probe_groups)
probe_groups$address <- as.numeric(probe_groups$address)

# Merge with your intensity data
ctrlData <- ctrlData %>%
  left_join(probe_groups, by = "address")

# Add probe group
ctrlData <- data.frame(ctrlData %>%
  mutate(  probe_group = case_when(
           color_label %in% green_channel_label & probe_type %in% paste0("C",1:6) ~ "Green Channel High",
           color_label %in% green_channel_label & probe_type %in% paste0("U",1:6) ~ "Green Channel Background",
           color_label %in% red_channel_label & probe_type %in% paste0("C",1:6) ~ "Red Channel High",
           color_label %in% red_channel_label & probe_type %in% paste0("U",1:6) ~ "Red Channel Background"
         )))

# Create the box plot to check the 
# Measure only in one channel, as specified before
df = ctrlData[ctrlData$channel == unlist(lapply(strsplit(ctrlData$probe_group, " "), function(x) x[1])),]

summarized_df <- df %>%
  select(slide, array, probe_group, value) %>%
  group_by(slide, array, probe_group) %>%
  summarise(mean_value = mean(value, na.rm = TRUE), .groups = 'drop')

sample_dependent_controls_bisulfite_conversion_I_array_vs_slide_scatter_plot <- ggplot(data = summarized_df, aes(x = as.factor(array), y = mean_value, color = probe_group)) +
  geom_point(size = 0.5) + # Made points very small
  labs(
    title = "Mean Value by Slide, Array, and Probe Group",
    y = "Mean Value",
    color = "Probe Group"
  ) +
  scale_color_manual(values = c(
    "Red Channel High" = "red",
    "Red Channel Background" = "lightpink",
    "Green Channel High" = "green",
    "Green Channel Background" = "lightgreen"
  )) +
  theme_minimal() +
  theme(
    axis.title.x = element_blank(), # Removed x-axis label
    panel.grid = element_blank(), # Removed all grid lines
    plot.title = element_text(size = 8), # Small plot title
    axis.text = element_text(size = 6), # Small axis text
    axis.title.y = element_text(size = 7), # Small y-axis title
    legend.title = element_text(size = 7), # Small legend title
    legend.text = element_text(size = 6), # Small legend text
    strip.text = element_text(size = 7) # Small facet labels
  ) +
  facet_wrap(~slide, scales = "free_x")


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

sample_dependent_controls_bisulfite_conversion_I_well <- 
  ggplot(df, aes(x = sample_well, y = value, fill = probe_group)) +
  geom_boxplot(position = position_dodge(0.8), width = 0.6) +
  facet_wrap(~sample_well, scales = "free_x", ncol = 8) + 
  scale_fill_manual(values = c(
    "Green Channel High" = "darkgreen",
    "Green Channel Background" = "lightgreen",
    "Red Channel High" = "darkred",
    "Red Channel Background" = "pink"
  )) +
  labs(
    title = "Bisulfite Conversion Controls type I by well",
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
sample_dependent_controls_bisulfite_conversion_I <- ggplot(ctrlData, aes(x = green, y = red, color = probe_type)) +
  # Add diagonal x=y line first (so it's behind points)
  geom_abline(intercept = 0, slope = 1, color = "gray50", linetype = "dashed") +
  # Plot points
  geom_point(size = 1, alpha = 0.3) +
  # Make axes equal and set limits
  coord_equal(xlim = c(5,15), ylim = c(5, 15)) +
  labs(
    title = "Bisulfite Conversion Controls type I",
    subtitle = "Each sample is represented by 12 points (C1,C2,C3,C4,C5,U1,U2,U3,U4,U5)",
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
  filter(probe_group == "Green Channel Background") %>%
  pull(green) %>%
  .[is.finite(.)] %>%
  mean(na.rm = TRUE)

# U4-U6 represent background for red channel
bkg_red <- ctrlData %>%
  filter(probe_group == "Red Channel Background") %>%
  pull(red) %>%
  .[is.finite(.)] %>%
  mean(na.rm = TRUE)

# Identify samples where controls (C1-C6) are below background thresholds
bad_samples <- ctrlData %>%
  group_by(sample) %>%
  summarize(
    green_low = mean(green[probe_group == "Green Channel High"], na.rm = TRUE) < bkg_green,
    red_low = mean(red[probe_group == "Red Channel High"], na.rm = TRUE) < bkg_red,
    .groups = "drop"
  ) %>%
  filter(green_low | red_low) %>%
  distinct(sample)

# Write output
write.csv(
  na.omit(bad_samples),
  file = opt$bs1_output,
  row.names = FALSE
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
ctrlData$sample_plate <- RGset@colData[ctrlData$sample,]$Plate
ctrlData$sample_well <- RGset@colData[ctrlData$sample,]$Well
ctrlData$slide <- RGset@colData[ctrlData$sample,]$Slide
ctrlData$array <- RGset@colData[ctrlData$sample,]$Array

# Comparing intensity distribution platewise
# Everything is evaluated in red channel so remove green
df = ctrlData[ctrlData$channel == "Red",]


summarized_df <- df %>%
  select(slide, array, value) %>%
  group_by(slide, array) %>%
  summarise(mean_value = mean(value, na.rm = TRUE), .groups = 'drop')


# Then, create the scatter plot using the summarized data with the requested changes
sample_dependent_controls_bisulfite_conversion_II_array_vs_slide_scatter_plot <- ggplot(data = summarized_df, aes(x = as.factor(array), y = mean_value)) +
  geom_point(size = 0.5) + # Made points very small
  labs(
    title = "Mean Value by Slide, Array, and Probe Group",
    y = "Mean Value"
  ) +
  theme_minimal() +
  theme(
    axis.title.x = element_blank(), # Removed x-axis label
    panel.grid = element_blank(), # Removed all grid lines
    plot.title = element_text(size = 8), # Small plot title
    axis.text = element_text(size = 6), # Small axis text
    axis.title.y = element_text(size = 7), # Small y-axis title
    legend.title = element_text(size = 7), # Small legend title
    legend.text = element_text(size = 6), # Small legend text
    strip.text = element_text(size = 7) # Small facet labels
  ) +
  facet_wrap(~slide, scales = "free_x")

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

  sample_dependent_controls_bisulfite_conversion_II_well <- 
  ggplot(df, aes(y = sample_well, x = value)) + 
  geom_boxplot(width = 0.6,  fill = "lightblue") +  # Single color
  labs(
    title = "Bisulfite Conversion Controls type II by well",
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
sample_dependent_controls_bisulfite_conversion_II <- ggplot(ctrlData, aes(x = red)) +
  geom_histogram(
    binwidth = 0.5,              
    fill = "gray80",              
    color = "white",              
    alpha = 0.8                    
  ) +
  labs(
    title = "Bisulfite Conversion Controls type II",
    subtitle = "Controls to be evaluated in the red channel\nEach sample is represented by 4 points",
    x = "Log2 Red Channel Intensity",
    y = "Count"
  ) +
  scale_x_continuous(
    limits = c(5, 15),           
    breaks = seq(5, 15, by = 1)    
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    panel.grid.major.x = element_line(color = "gray90"),
    panel.grid.minor.x = element_blank(),
    aspect.ratio = 0.6            
  )
  
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
  row.names = FALSE
)

# Finally we will check if there is any plate specifc bias in this step. It coudl be that all the samples we have flagged are actually coming from the same plate
pdf(opt$pdf_output, width = 8, height = 11)  # Increased height
grid.draw(ggplotGrob(sample_dependent_controls_bisulfite_conversion_I))
grid.newpage()
grid.draw(ggplotGrob(sample_dependent_controls_bisulfite_conversion_I_plate))
grid.newpage()
grid.draw(ggplotGrob(sample_dependent_controls_bisulfite_conversion_I_well))
grid.newpage()
grid.draw(ggplotGrob(sample_dependent_controls_bisulfite_conversion_I_array_vs_slide_scatter_plot))
grid.newpage()
grid.draw(ggplotGrob(sample_dependent_controls_bisulfite_conversion_II))
grid.newpage()
grid.draw(ggplotGrob(sample_dependent_controls_bisulfite_conversion_II_plate))
grid.newpage()
grid.draw(ggplotGrob(sample_dependent_controls_bisulfite_conversion_II_well))
grid.newpage()
grid.draw(ggplotGrob(sample_dependent_controls_bisulfite_conversion_II_array_vs_slide_scatter_plot))
dev.off()


