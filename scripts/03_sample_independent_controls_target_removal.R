#!/usr/bin/env Rscript
# =============================================================================
# 03_sample_independent_controls_target_removal.R
# Author: Vartika Bisht
# Created: 13-Aug-2025
# Description:
# This script analyzes target removal control probes to evaluate the efficiency
# of the post-extension stripping step. These controls are designed to be
# washed off the array, and therefore, a low signal intensity is expected,
# indicating successful removal of non-specifically bound targets.
# Arguments:
#    --rgset                → Path to the RGChannelSet object file (.RData)
#    --platform             → DNA methylation platform (e.g., EPIC, 450K)
#    --manifestkey          → Path to the manifest key file (default: data/manifest.annotation.key.csv)
#    --assembly             → Genome assembly version
#    --flagged              → Output CSV file for flagged target removal failures
#    --pdf_output           → Output PDF file for plots
# Usage:
#   Rscript 03_sample_independent_controls_target_removal.R \
#     --rgset <FILE.RData> --platform <PLATFORM> --assembly <ASSEMBLY> \
#     --flagged <target_removal_failures.csv> --pdf_output <output.pdf>
# Notes:
#   - Target removal is a critical wash step that removes weakly bound or unbound
#     DNA, helping to reduce background noise.
#   - The two target removal control probes are expected to have low signal
#     intensity, indicating that they were efficiently washed off the array.
#   - Channel-specific Evaluation: The performance of these controls is
#     monitored exclusively in the green channel.
#   - Flagging Condition: As there is no standard definition for a "low" signal,
#     samples are flagged if the mean intensity of these two probes is more than two
#     standard deviations above the median intensity of all samples.
# =============================================================================

# Load R packages silently
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
              help="Output file for target removal failures"),
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

# | Purpose        | Name           | Number on the Array | Evaluate Green (GRN) | Evaluate Red (RED) | Expected Intensity |
# | :------------- | :------------- | :------------------ | :------------------- | :----------------- | :----------------- |
# | Target Removal | Target Removal | 1, 2                | +                    | -                  | Low                |

ctrl = "TARGET REMOVAL"
red <- getRed(RGset)
green <- getGreen(RGset)
ctrlAddress <- getControlAddress(RGset, controlType = ctrl)
ctrlData <-  data.frame(cbind(channel = "Green", melt(green[ctrlAddress, ], varnames = c("address", "sample"))))
ctrlData$value <- log2(ctrlData$value)
ctrlData$slide <- RGset@colData[ctrlData$sample,]$Slide
ctrlData$array <- RGset@colData[ctrlData$sample,]$Array


probe_groups <- controlprobes[controlprobes$Type == ctrl,]
colnames(probe_groups) = c("address","control_probe","hypothetical_color","name")
probe_groups <- data.frame(probe_groups)
probe_groups$address <- as.numeric(probe_groups$address)

# Merge with your intensity data
ctrlData <- ctrlData %>%
  left_join(probe_groups, by = "address")


# Define bad samples by mean + 2 * sd
bkg <-  mean(ctrlData$value, na.rm = TRUE) + 2 * sd(ctrlData$value, na.rm = TRUE)
bad_samples <- ctrlData[ctrlData$value > bkg, ]$sample

# Write output
write.csv(
  na.omit(bad_samples),
  file = opt$flagged,
  row.names = FALSE
)


sample_independent_controls_target_removal_by_array <- ggplot(ctrlData, aes(x = value, y = array), color = "darkgreen") + # Map color here for box outlines
  geom_boxplot(show.legend = FALSE) + # Use geom_boxplot and hide the legend for color
  scale_x_continuous(limits = c(5, 17)) +
  labs(
    title = "Control Probe Intensities by Array (Box Plots)", # Updated plot title
    x = "Intensity Value",                      # X-axis label
    y = "Array ID"                              # Y-axis label
  ) +
  theme_minimal() +                               # Use a clean, minimal theme
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold") # Center and bold the plot title
  )


summarized_df <- data.frame(ctrlData %>%
  group_by(slide, name) %>%
  summarise(value = median(value), .groups = 'drop'))  %>%
  mutate(slide = as.factor(slide))
summarized_df <- summarized_df[is.finite(summarized_df$value),]

sample_independent_controls_target_removal_by_slide <- ggplot(summarized_df, aes(y = slide, x = value) , color = "darkgreen") + # Swap x and y aesthetics
  geom_point(size = 3) + 
  labs(
    title = "Target removal controls by slide",
    x = "Median Intensity", # Update x-axis label
    y = "Slide ID"    # Update y-axis label
  ) +
  theme_minimal() + # A clean, minimal theme
  scale_x_continuous(limits = c(5, 17)) +  
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold")
  )


sample_independent_controls_target_removal_cutoff <- ggplot(ctrlData, aes(x = name, y = value)) +
  geom_boxplot(
    outlier.color = "black",  
    alpha = 0.8 
  ) +
  geom_hline(yintercept = bkg, linetype = "dashed", color = "blue") +  
  labs(
    title = "Target Removal Control",
    x = "",
    y = "Log2 Intensity in green channel",
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5),
    legend.position = "bottom",
    panel.grid.major.x = element_blank()  # Remove vertical grid lines for cleaner look
  ) 


pdf(opt$pdf_output, width = 8, height = 11)
grid.draw(ggplotGrob(sample_independent_controls_target_removal_by_array))
grid.newpage()
grid.draw(ggplotGrob(sample_independent_controls_target_removal_by_slide))
grid.newpage()
grid.draw(ggplotGrob(sample_independent_controls_target_removal_cutoff))
dev.off()
