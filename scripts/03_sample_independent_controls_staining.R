#!/usr/bin/env Rscript
# =============================================================================
# 03_sample_independent_controls_staining.R
# Author: Vartika Bisht
# Created: 13-Aug-2025
# Description:
# This script analyzes staining controls to assess the efficiency of the
# staining process in both the red and green channels. These controls are
# independent of the hybridization and extension steps, providing a direct
# measure of staining quality and reagent performance.
# Arguments:
#    --rgset                → Path to the RGChannelSet object file (.RData)
#    --platform             → DNA methylation platform (e.g., EPIC, 450K)
#    --manifestkey          → Path to the manifest key file (default: data/manifest.annotation.key.csv)
#    --assembly             → Genome assembly version
#    --flagged              → Output CSV file for flagged staining failures
#    --pdf_output           → Output PDF file for plots
# Usage:
#   Rscript 03_sample_independent_controls_staining.R \
#     --rgset <FILE.RData> --platform <PLATFORM> --assembly <ASSEMBLY> \
#     --flagged <staining_failures.csv> --pdf_output <output.pdf>
# Notes:
#   - Staining controls are sample-independent and are designed to evaluate the
#     effectiveness of the staining reagents used in the assay.
#   - The controls consist of two pairs of probes:
#     - DNP-immobilized beads are assessed in the red channel.
#     - Biotin-immobilized beads are assessed in the green channel.
#   - For both channels, there are high-intensity probes and background probes.
#   - Expected Performance: For each channel, the signal from the high-intensity
#     probe is expected to be significantly greater than that of the background probe
#     (i.e., High > Background).
#   - Flagging Condition: The script flags any sample where this expected relationship
#     is not met, as it may indicate an issue with the staining process.
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
              help="Genome assembly version", metavar="character"),
  make_option(c("-f", "--flagged"), type="character", default=NULL,
              help="Output file for staining failures"),
  make_option(c("-o", "--pdf_output"), type="character", default=NULL,
              help="Output PDF file for plots")
)


# Parse command line options
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)


# VALIDATE REQUIRED ARGUMENTS
required_args <- c("rgset", "platform", "manifestkey", "assembly" , "flagged" , "pdf_output")

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

# | Staining Type  | Replicate | Evaluate Green (GRN) | Evaluate RED (RED) | Signal Level  |
# | :------------- | :-------- | :------------------- | :----------------- | :------------ |
# | DNP (High)     | 1         | -                    | +                  | High          |
# | DNP (Bgnd)     | 1         | -                    | +                  | Background    |
# | Biotin (Med)   | 1         | +                    | -                  | High          |
# | Biotin (Bgnd)  | 1         | +                    | -                  | Background    |

# --- Sample Independent Controls (STAINING) ---
# These are specific probes used for staining quality assessment:
#   - 27630314, STAINING, Red, DNP (High),
#   - 43603326, STAINING, Purple, DNP (Bkg),
#   - 41666334, STAINING, Green, Biotin (High),
#   - 34648333, STAINING, Blue, Biotin (Bkg),


ctrl = "STAINING"
red <- getRed(RGset)
green <- getGreen(RGset)
ctrlAddress <- getControlAddress(RGset, controlType = ctrl)
ctrlData <-  rbind(
            cbind(channel = "Red", melt(red[ctrlAddress, ], varnames = c("address", "sample"))),
            cbind(channel = "Green", melt(green[ctrlAddress, ], varnames = c("address", "sample"))))
ctrlData$value <- log2(ctrlData$value)
ctrlData$slide <- RGset@colData[ctrlData$sample,]$Slide
ctrlData$array <- RGset@colData[ctrlData$sample,]$Array

rm(red, green, ctrlAddress)

probe_groups <- controlprobes[controlprobes$Type == ctrl,]
probe_groups <- probe_groups[probe_groups$Color != -99,]
probe_groups$expected_channel <- ifelse(grepl("DNP", probe_groups$ExtendedType), "Red", "Green")
colnames(probe_groups) = c("address","control_probe","hypothetical_color","name","expected_channel")
probe_groups <- data.frame(probe_groups)
probe_groups$address <- as.numeric(probe_groups$address)

# Merge with your intensity data
ctrlData <- ctrlData %>%
  left_join(probe_groups, by = "address")

# Pivot to wide format (one row per address-sample pair)
ctrlDatal <- ctrlData %>%
  pivot_wider(
    names_from = channel,      # Column to split (Red/Green)
    values_from = value,       # Values to fill into new columns
    names_glue = "{tolower(channel)}"  # Force lowercase column names
  )

ctrlDatal <- data.frame(ctrlDatal)

# Probe Type Biotin (Bkg) Biotin (High) DNP (Bkg) DNP (High)
sample_independent_controls_staining <- ggplot(ctrlDatal, aes(x = green, y = red, color = name)) +
  # Add diagonal x=y line first (so it's behind points)
  geom_abline(intercept = 0, slope = 1, color = "gray50", linetype = "dashed") +
  # Plot points
  geom_point(size = 1, alpha = 0.3) +
  # Make axes equal and set limits
  coord_equal(xlim = c(5,15), ylim = c(5, 15)) +
  labs(
    title = "Staining controls",
    x = "Log2 Green Channel Intensity",
    y = "Log2 Red Channel Intensity",
    color = "Probe Type"
  ) +
  theme_minimal() + 
  scale_color_manual(values = c(
    "DNP (High)" = "darkgreen",
    "DNP (Bkg)" = "lightgreen",
    "Biotin (High)" = "darkred",
    "Biotin (Bkg)" = "pink"
  )) +
  theme(
    legend.position = "bottom",
    legend.title = element_text(size = 10),
    plot.title = element_text(hjust = 0.5),
    aspect.ratio = 1  # Ensure square plot area
  ) +
  guides(color = guide_legend(override.aes = list(size = 3, alpha = 1)))

rm(ctrlDatal)
gc()

# Keep probe intensity value enteries from the expected channels
ctrlData <- ctrlData[ctrlData$channel == ctrlData$expected_channel,]

sample_independent_controls_staining_by_array <- ggplot(ctrlData, aes(x = array, y = value, fill = name)) +
  geom_boxplot(position = position_dodge(0.8), width = 0.6) +
  facet_wrap(~array, scales = "free_x", ncol = 4) + 
  scale_fill_manual(values = c(
    "DNP (High)" = "darkgreen",
    "DNP (Bkg)" = "lightgreen",
    "Biotin (High)" = "darkred",
    "Biotin (Bkg)" = "pink"
  )) +
  labs(
    title = "Staining controls by position on slide (aka array)",
    x = "Position on slide",
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

sample_independent_controls_staining_by_slide <- ggplot(ctrlData,  aes(x = name, y = value, fill = name ) ) +
  geom_boxplot() + 
  facet_wrap(~ slide, scales = "free_x") +
  scale_fill_manual( 
    values = c("Biotin" = "darkgreen", "DNP" = "darkred"), 
    name = "Probe Type" ) +
  labs(
    title = "Staining Controls by Slide (Box Plots)",
    x = "Probe Type", 
    y = "Intensity" ) +
  theme_minimal() + 
  scale_fill_manual(values = c(
    "DNP (High)" = "darkgreen",
    "DNP (Bkg)" = "lightgreen",
    "Biotin (High)" = "darkred",
    "Biotin (Bkg)" = "pink"
  )) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1) 
  )

bad_samples = data.frame("V1"=union(ctrlData[ctrlData$name == "DNP (High)" & ctrlData$value < median(ctrlData[ctrlData$name == "DNP (Bkg)",]$value),]$sample,ctrlData[ctrlData$name == "Biotin (High)" & ctrlData$value < median(ctrlData[ctrlData$name == "Biotin (Bkg)",]$value),]$sample))


# Write output
write.csv(
  na.omit(bad_samples),
  file = opt$flagged,
  row.names = FALSE
)


pdf(opt$pdf_output, width = 8, height = 11)
grid.draw(ggplotGrob(sample_independent_controls_staining))
grid.newpage()
grid.draw(ggplotGrob(sample_independent_controls_staining_by_array))
grid.newpage()
grid.draw(ggplotGrob(sample_independent_controls_staining_by_slide))
dev.off()


