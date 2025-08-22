#!/usr/bin/env Rscript
# =============================================================================
# 03_sample_dependent_controls_specificity.R
# Author: Vartika Bisht
# Created: 13-Aug-2025
# Description:
# This script analyzes specificity controls to monitor non-specific primer
# extension and proper base incorporation for both Infinium I and II probes.
# The analysis evaluates probe performance in both the red and green channels
# to identify samples with potential assay failures.
# Arguments:
#    --rgset                → Path to the RGChannelSet object file (.RData)
#    --platform             → DNA methylation platform (e.g., EPIC, 450K)
#    --manifestkey          → Path to the manifest key file (default: data/manifest.annotation.key.csv)
#    --assembly             → Assembly version of the DNA methylation platform
#    --sp1_output           → Output CSV file path for SPECIFICITY I failures
#    --sp2_output           → Output CSV file path for SPECIFICITY II failures
#    --pdf_output           → Output PDF file path for QC plots
# Usage:
#   Rscript 03_sample_dependent_controls_specificity.R \
#     --rgset <FILE.RData> --platform <PLATFORM> --assembly <ASSEMBLY> \
#     --sp1_output <sp1_failures.csv> --sp2_output <sp2_failures.csv> \
#     --pdf_output <output.pdf>
# Notes:
#   - Specificity controls are designed to monitor primer annealing and base incorporation
#     specificity, identifying non-specific binding or extension.
#   - Specificity I Controls (Infinium I Probes):
#      - These probes come in Perfect Match (PM) and Mismatch (MM) pairs.
#      - PM probes (A/T match) should have high signal, while MM probes (G/T mismatch)
#         should have low signal.
#      - The analysis evaluates PM/MM pairs in both the red and green channels.
#      - Flagging Condition: Samples are flagged if the PM signal is lower than the
#         MM signal in either channel.
#   - Specificity II Controls (Infinium II Probes):
#      - There are 3 controls designed to monitor extension specificity.
#      - All 3 probes are measured in the red channel and are expected to show
#         high intensity.
#      - Flagging Condition: Samples are flagged if the signal intensity of any
#         of these probes falls below an empirical cutoff of mean intensity - (2 * standard deviation).
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
  make_option(c("-o", "--sp1_output"), type="character", default=NULL,
              help="Output file for SPECIFICITY I failures"),
  make_option(c("-t", "--sp2_output"), type="character", default=NULL,
              help="Output file for SPECIFICITY II failures"),
  make_option(c("-d", "--pdf_output"), type="character", default=NULL,
              help="Output PDF file for plots")
)

# Parse command line options
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)


# VALIDATE REQUIRED ARGUMENTS
required_args <- c("rgset", "platform", "manifestkey", "assembly" , "sp1_output" , "sp2_output", "pdf_output")

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

# Load the RGset data
message("Loading RGset data from: ", opt$rgset)
load(opt$rgset)

# ==================== SPECIFICITY I Controls ====================
# Control Probe Details:
# ======================
# * Specificity I Controls (12 probes - measured in both Green and Red channels):
#   - PM (Perfect Match) probes end with A (A/T match - should have HIGH signal)
#   - MM (Mismatch) probes end with G (G/T mismatch - should have LOW signal)
#   - PM1-3/MM1-3 measured in Red channel
#   - PM4-6/MM4-6 measured in Green channel

ctrl = "SPECIFICITY I"
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
rm(red, green, ctrlAddress)

# Perfect match (PM) and mismatch (MM) probes
# PM should be higher than MM
probe_groups <- controlprobes[controlprobes$Type == ctrl,]
probe_groups$group <- unlist(lapply(strsplit(probe_groups$ExtendedType, " "), function(x) x[3]))
probe_groups$ExtendedType <- substr(probe_groups$ExtendedType,16,17)
probe_groups$match <- paste(probe_groups$ExtendedType,probe_groups$group, sep = "")
colnames(probe_groups) <- c("address", "control_type", "color_label", "match_status", "group", "match")
probe_groups <- data.frame(probe_groups)
probe_groups$address <- as.numeric(probe_groups$address)

# Merge with your intensity data
ctrlData <- ctrlData %>%
  left_join(probe_groups, by = "address")

# Expected channels for each match type:
# Add probe group
ctrlData <- data.frame(ctrlData %>%
  mutate(  probe_group = case_when(
           color_label %in% green_channel_label & match %in% paste0("PM",1:6) ~ "Green Channel High",
           color_label %in% green_channel_label & match %in% paste0("MM",1:6) ~ "Green Channel Background",
           color_label %in% red_channel_label & match %in% paste0("PM",1:6) ~ "Red Channel High",
           color_label %in% red_channel_label & match %in% paste0("MM",1:6) ~ "Red Channel Background"
         )))


# Create the box plot to check the 
# Measure only in one channel, as specified before
df = ctrlData[ctrlData$channel == unlist(lapply(strsplit(ctrlData$probe_group, " "), function(x) x[1])),]


summarized_df <- df %>%
  select(slide, array, probe_group, value) %>%
  group_by(slide, array, probe_group) %>%
  summarise(mean_value = mean(value, na.rm = TRUE), .groups = 'drop')


sample_dependent_controls_specificity_I_array_vs_slide_scatter_plot <- ggplot(data = summarized_df, aes(x = as.factor(array), y = mean_value, color = probe_group)) +
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


sample_dependent_controls_specificity_I_plate <- 
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
    title = "Specificity Controls type I by plate",
    x = "Sample Plate",
    y = "Log2 Intensity",
    fill = "Match Type"
  ) +
  theme_minimal() +
  theme(
    axis.title.x = element_blank(),     
    axis.text.x = element_blank(),      
    axis.ticks.x = element_blank(),   
    strip.background = element_blank(),  
    panel.spacing = unit(0.5, "lines")  
  )

sample_dependent_controls_specificity_I_well <- 
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
    title = "Specificity Controls type I by well",
    x = "Sample Plate",
    y = "Log2 Intensity",
    fill = "Match Type"
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
  row.names = FALSE
)

# ==================== SPECIFICITY II Controls ====================

# * Specificity II Controls (3 probes - measured in Red channel only):
#   - Should incorporate "A" base across nonpolymorphic T (Red channel signal)
#   - Non-specific "G" base incorporation reduces Red channel signal

ctrl = "SPECIFICITY II"
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
sample_dependent_controls_specificity_II_array_vs_slide_scatter_plot <- ggplot(data = summarized_df, aes(x = as.factor(array), y = mean_value)) +
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

sample_dependent_controls_specificity_II_plate <- 
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

  sample_dependent_controls_specificity_II_well <- 
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
  row.names = FALSE
)

# Finally we will check if there is any plate specifc bias in this step. It coudl be that all the samples we have flagged are actually coming from the same plate
pdf(opt$pdf_output, width = 8, height = 11) 
grid.draw(ggplotGrob(sample_dependent_controls_specificity_I))
grid.newpage()
grid.draw(ggplotGrob(sample_dependent_controls_specificity_I_plate))
grid.newpage()
grid.draw(ggplotGrob(sample_dependent_controls_specificity_I_well))
grid.newpage()
grid.draw(ggplotGrob(sample_dependent_controls_specificity_I_array_vs_slide_scatter_plot))
grid.newpage()
grid.draw(ggplotGrob(sample_dependent_controls_specificity_II))
grid.newpage()
grid.draw(ggplotGrob(sample_dependent_controls_specificity_II_plate))
grid.newpage()
grid.draw(ggplotGrob(sample_dependent_controls_specificity_II_well))
grid.newpage()
grid.draw(ggplotGrob(sample_dependent_controls_specificity_II_array_vs_slide_scatter_plot))
dev.off()