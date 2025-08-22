#!/usr/bin/env Rscript
# =============================================================================
# 03_sample_independent_controls_hybridization.R
# Author: Vartika Bisht
# Created: 13-Aug-2025
# Description:
# This script analyzes hybridization control probes to assess the overall
# performance of the methylation assay. These controls use synthetic targets
# at varying concentrations to verify the efficiency of the hybridization
# process itself, providing a sample-independent quality check.
# Arguments:
#    --rgset                → Path to the RGChannelSet object file (.RData)
#    --platform             → DNA methylation platform (e.g., EPIC, 450K)
#    --manifestkey          → Path to the manifest key file (default: data/manifest.annotation.key.csv)
#    --assembly             → Genome assembly version
#    --flagged              → Output CSV file for flagged hybridization failures
#    --pdf_output           → Output PDF file for plots
# Usage:
#   Rscript 03_sample_independent_controls_hybridization.R \
#     --rgset <FILE.RData> --platform <PLATFORM> --assembly <ASSEMBLY> \
#     --flagged <hyb_failures.csv> --pdf_output <output.pdf>
# Notes:
#   - Hybridization controls use synthetic targets at low, medium, and high concentrations
#     to test the efficiency of the hybridization process. The controls are sample-independent.
#   - Expected Performance: The probe intensities are expected to follow a linear
#     response corresponding to the concentration gradient: High > Medium > Low.
#   - Channel-specific Evaluation: The performance of these controls is monitored
#     exclusively in the green channel.
#   - Flagging Condition: A sample is flagged if the intensity order is not as expected
#     (i.e., if High < Medium or Medium < Low). This indicates a potential failure
#     in the hybridization step of the assay.
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
              help="Output file for hybridization failures"),
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

# --- Hybridization Controls Design --- 450K
# These controls specifically assess the efficiency of the hybridization process itself.

# | Purpose       | Name        | Number on the Array | Evaluate Green (GRN) | Evaluate Red (RED) | Expected Intensity |
# | :------------ | :---------- | :------------------ | :------------------- | :----------------- | :----------------- |
# | Hybridization | Hyb (Low)   | 1                   | +                    | -                  | Low                |
# | Hybridization | Hyb (Medium)| 1                   | +                    | -                  | Medium             |
# | Hybridization | Hyb (High)  | 1                   | +                    | -                  | High               |

# --- Specific Hybridization Control Probes ---
# These are the specific probe IDs utilized for assessing hybridization quality:
#   - 28684356, HYBRIDIZATION, Green, Hyb (High)
#   - 26772442, HYBRIDIZATION, Blue, Hyb (Medium)
#   - 21771417, HYBRIDIZATION, Black, Hyb (Low)


ctrl = "HYBRIDIZATION"
red <- getRed(RGset)
green <- getGreen(RGset)
ctrlAddress <- getControlAddress(RGset, controlType = ctrl)
# Everything is to be measured in the green channel
ctrlData <-  data.frame(cbind(channel = "Green", melt(green[ctrlAddress, ], varnames = c("address", "sample"))))
ctrlData$value <- log2(ctrlData$value)
ctrlData$slide <- RGset@colData[ctrlData$sample,]$Slide
ctrlData$array <- RGset@colData[ctrlData$sample,]$Array

probe_groups <- controlprobes[controlprobes$Type == ctrl,]
probe_groups$ExtendedType <- substr(probe_groups$ExtendedType,6,nchar(probe_groups$ExtendedType)-1)
colnames(probe_groups) = c("address","control_probe","hypothetical_color","name")
probe_groups <- data.frame(probe_groups)
probe_groups$address <- as.numeric(probe_groups$address)

# Merge with your intensity data
ctrlData <- ctrlData %>%
  left_join(probe_groups, by = "address")

# Everything is measured in the green channel and should be low
# Pivot to wide format (one row per address-sample pair)
ctrlDatal <- ctrlData %>%
  pivot_wider(
    names_from = channel,      # Column to split (Red/Green)
    values_from = value,       # Values to fill into new columns
    names_glue = "{tolower(channel)}"  # Force lowercase column names
  )

ctrlDatal <- data.frame(ctrlDatal)

sample_independent_controls_hybridization <- ggplot(ctrlDatal, aes(x = green, fill = name)) +
  # Plot points
  geom_density() +
  scale_fill_manual(
    values = c(
      "High" = "#006400",    # DarkGreen
      "Medium" = "#32CD32",  # LimeGreen
      "Low" = "#90EE90"      # LightGreen
    )
  ) +
  labs(
    title = "Hybridization controls",
    x = "Log2 Green Channel Intensity",
    y = "Denisty",
    fill = "Intensity Group"
  ) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    legend.title = element_text(size = 10),
    plot.title = element_text(hjust = 0.5),
    aspect.ratio = 1  # Ensure square plot area
  ) +
  guides(color = guide_legend(override.aes = list(size = 3, alpha = 1))) + coord_equal(xlim = c(10,16))

rm(ctrlDatal)
gc()

ctrlData$name <- factor(ctrlData$name, levels = c("Low", "Medium", "High"))
ctrlData$array <- as.factor(ctrlData$array)

sample_independent_controls_hybridization_by_array <- ggplot(ctrlData, aes(x = array, y = value, fill = name)) +
  geom_boxplot(position = position_dodge(0.8), width = 0.6) +
  facet_wrap(~array, scales = "free_x", ncol = 4) +
  scale_fill_manual(
    values = c(
      "High" = "#006400",    # DarkGreen
      "Medium" = "#32CD32",  # LimeGreen
      "Low" = "#90EE90"      # LightGreen
    ),
    # Add 'breaks' to explicitly set the order in the legend,
    # though factor levels usually handle this automatically.
    breaks = c("Low", "Medium", "High"),
    name = "Probe Type" # Keep the legend title
  ) +
  labs(
    title = "Hybridization controls by position on slide (aka array)",
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

sample_independent_controls_hybridization_by_slide <- ggplot(
  ctrlData, # Use the original (filtered) data to show distributions and outliers
  aes(x = name, y = value, fill = name) # 'name' on x-axis, 'value' on y-axis, fill by 'name'
) +
  geom_boxplot() + # Use geom_boxplot for box plots
  facet_wrap(~ slide, scales = "free_x") + # Facet by 'slide'
  scale_fill_manual( # Use scale_fill_manual for box plot fill
    values = c(
      "High" = "#006400",  # DarkGreen
      "Medium" = "#32CD32", # LimeGreen
      "Low" = "#90EE90"    # LightGreen
    ),
    breaks = c("Low", "Medium", "High"), # Ensure legend order
    name = "Probe Type" # Legend title
  ) +
  labs(
    title = "Target Removal Controls by Slide (Box Plots)",
    x = "Probe Type", # Update x-axis label
    y = "Intensity" # Update y-axis label
  ) +
  theme_minimal() + # A clean, minimal theme
  scale_y_continuous(limits = c(5, 17)) + # Set y-axis limits (corresponding to original x-axis limits)
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"), # Center and bold the title
    axis.text.x = element_text(angle = 45, hjust = 1) # Rotate x-axis labels for readability
  )



# Define bad samples by the samples which do not follow high medium low pattenr
wide_ctrl_data <- ctrlData %>%
  select(sample, value, name) %>%
  pivot_wider(names_from = name, values_from = value)
bad_samples <- wide_ctrl_data %>%
  filter(!(High > Medium & Medium > Low)) %>% # Filter for rows where the condition is FALSE
  select(sample) %>% # Select only the 'sample' column
  distinct() # In case a sample could be flagged by multiple conditions, ensure unique names


# Write output
write.csv(
  data.frame(bad_samples),
  file = opt$flagged,
  row.names = FALSE
)


pdf(opt$pdf_output, width = 8, height = 11)
grid.draw(ggplotGrob(sample_independent_controls_hybridization))
grid.newpage()
grid.draw(ggplotGrob(sample_independent_controls_hybridization_by_array))
grid.newpage()
grid.draw(ggplotGrob(sample_independent_controls_hybridization_by_slide))
dev.off()
