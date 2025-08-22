#!/usr/bin/env Rscript
# =============================================================================
# 03_sample_independent_controls_extension.R
# Author: Vartika Bisht
# Created: 13-Aug-2025
# Description:
# This script analyzes the extension control probes to assess the efficiency
# of the single-base extension reaction for each nucleotide (A, T, C, G).
# As these controls use a hairpin probe design, they are sample-independent
# and provide a crucial check on the overall reagent and BeadChip performance.
# Arguments:
#    --rgset                → Path to the RGChannelSet object file (.RData)
#    --platform             → DNA methylation platform (e.g., EPIC, 450K)
#    --manifestkey          → Path to the manifest key file (default: data/manifest.annotation.key.csv)
#    --assembly             → Genome assembly version
#    --flagged              → Output CSV file for flagged extension failures
#    --pdf_output           → Output PDF file for plots
# Usage:
#   Rscript 03_sample_independent_controls_extension.R \
#     --rgset <FILE.RData> --platform <PLATFORM> --assembly <ASSEMBLY> \
#     --flagged <extension_failures.csv> --pdf_output <output.pdf>
# Notes:
#   - Extension controls are "sample-independent" because they use a hairpin probe
#     design that contains its own complementary strand, eliminating the need for
#     genomic DNA as a template.
#   - These controls test the efficiency of the single-base extension of A, T, C, and G.
#   - Expected Performance:
#       - A and T probes are expected to have high intensity in the red channel.
#       - C and G probes are expected to have high intensity in the green channel.
#   - Flagging Condition: As there is no standard definition for "high" intensity,
#     a sample is flagged if the intensity of any extension probe is abnormally low.
#     The specific filter is defined as signal intensity < (mean intensity - 2 * standard deviation).
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
              help="Output file for extension failures"),
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

# Check if RGset exists after loading
if (!exists("RGset")) {
  stop("RGset object not found in the RData file. Please check the input file.")
}

# Expectation -- set up of illumina probes
# Extension Pair	Replicate	Evaluate Green (GRN)	Evaluate RED (RED)	Signal Level
# Extension (A), (T)	2	-	+	High
# Extension (C), (G)	2	+	-	High
# EXTENSION
# 63642461,EXTENSION,Red,Extension (A),
# 47640365,EXTENSION,Purple,Extension (T),
# 74666473,EXTENSION,Green,Extension (C),
# 31698466,EXTENSION,Blue,Extension (G),
# Extension control probe design

ctrl = "EXTENSION"
red <- getRed(RGset)
green <- getGreen(RGset)
ctrlAddress <- getControlAddress(RGset, controlType = ctrl)
ctrlData <-  rbind(
            cbind(channel = "Red", melt(red[ctrlAddress, ], varnames = c("address", "sample"))),
            cbind(channel = "Green", melt(green[ctrlAddress, ], varnames = c("address", "sample"))))
ctrlData$value <- log2(ctrlData$value)
ctrlData$slide <- RGset@colData[ctrlData$sample,]$Slide
ctrlData$array <- RGset@colData[ctrlData$sample,]$Array


probe_groups <- controlprobes[controlprobes$Type == ctrl,]
probe_groups$ExtendedType <- substr(probe_groups$ExtendedType,12,12)
probe_groups$expected_channel <- ifelse(probe_groups$ExtendedType %in% c("A","T"),"Red","Green")
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


# Get axis limits based on data range
sample_independent_controls_extension <- ggplot(ctrlDatal, aes(x = green, y = red, color = name)) +
  # Add diagonal x=y line first (so it's behind points)
  geom_abline(intercept = 0, slope = 1, color = "gray50", linetype = "dashed") +
  # Plot points
  geom_point(size = 1, alpha = 0.3) +
  # Make axes equal and set limits
  coord_equal(xlim = c(5,15), ylim = c(5, 15)) +
  labs(
    title = "Extension controls",
    x = "Log2 Green Channel Intensity",
    y = "Log2 Red Channel Intensity",
    color = "Probe Type"
  ) +
  theme_minimal() +
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

sample_independent_controls_extension_by_array <- ggplot(ctrlData, aes(x = array, y = value, fill = name)) +
  geom_boxplot(position = position_dodge(0.8), width = 0.6) +
  facet_wrap(~array, scales = "free_x", ncol = 4) +
  labs(
    title = "Extension controls by position on slide (aka array)",
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
  ) + scale_y_continuous(limits = c(6,15))

# Create the box plot
sample_independent_controls_extension_by_slide <- ggplot(
  ctrlData,
  aes(x = name, y = value, fill = expected_channel)
) +
  geom_boxplot() + 
  facet_wrap(~ slide, scales = "free_x") + 
  scale_fill_manual(
    values = c("Green" = "darkgreen", "Red" = "darkred"),
    name = "Channel" 
  ) +
  labs(
    title = "Extension Controls by Slide (Box Plots)",
    x = "Base",
    y = "Intensity" 
  ) +
  theme_minimal() + 
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1) 
  ) + scale_y_continuous(limits = c(6,15))



# Define bad samples by mean - 2 * sd
finite_values <- ctrlData$value[is.finite(ctrlData$value)]
mean_of_finite <- mean(finite_values, na.rm = TRUE)
sd_of_finite <- sd(finite_values, na.rm = TRUE)
bad_samples <- ctrlData[ctrlData$value < (mean_of_finite - 2 * sd_of_finite), ]$sample

# Write output
write.csv(
  unique(na.omit(bad_samples)),
  file = opt$flagged,
  row.names = FALSE
)


pdf(opt$pdf_output, width = 8, height = 11)
grid.draw(ggplotGrob(sample_independent_controls_extension))
grid.newpage()
grid.draw(ggplotGrob(sample_independent_controls_extension_by_array))
grid.newpage()
grid.draw(ggplotGrob(sample_independent_controls_extension_by_slide))
dev.off()


