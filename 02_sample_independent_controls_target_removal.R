# Rscript: 02_sample_independent_controls_target_removal.R
# Author: Vartika Bisht
# Date: 22 May 2025

# --- Purpose ---
# This script is specifically designed to evaluate the efficiency of the target removal (stripping) step
# that occurs after the extension reaction in the experimental workflow.

# --- Target Removal Controls Design ---
# These controls are engineered to test the effectiveness of the post-extension stripping step.
# During this process, control oligonucleotides are extended using the probe sequence as a template,
# which generates labeled targets. The probe sequences themselves are specifically designed
# to prevent unintended extension from the probe.

# | Purpose        | Name           | Number on the Array | Evaluate Green (GRN) | Evaluate Red (RED) | Expected Intensity |
# | :------------- | :------------- | :------------------ | :------------------- | :----------------- | :----------------- |
# | Target Removal | Target Removal | 1, 2                | +                    | -                  | Low                |

# --- Rationale for Analysis ---
# A crucial expectation for all target removal controls is that they should result in a low signal
# when compared to the hybridization controls. This low signal serves as a key indicator that
# the targets were efficiently removed following the extension reaction. These specific target
# removal controls are incorporated into the Hybridization Buffer RA1. It is important to note
# that the performance of these target removal controls should only be monitored in the green channel.

# --- Input Parameters ---
# -r / --RGset_file: Path to the RData file containing the RGset object. 
# -c / --controls_probes_file: Path to a CSV file that details the control probe IDs

# --- Output Files ---
# -p / --output_pdf: Path and filename for the generated PDF report. 
# -f / --flagged_csv: Path and filename for a CSV file that will list any flagged
#                     samples or arrays.

# --- Technical Notes ---
# This script includes provisions for silent library loading and robust error handling.

# --- Usage ---
# Rscript 02_sample_independent_controls_target_removal.R \
#   -r RGset.RData \
#   -c controls_probes.csv \
#   -p sample_independent_controls_target_removal_overview.pdf \
#   -f Flagged_target_removal.csv

# Load R packages silently
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
  make_option(c("-f", "--flag_output"), type="character", default=NULL,
              help="Output file for target removal failures"),
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


ctrl = "TARGET REMOVAL"
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

probe_groups <- control_probes[control_probes$V2 == ctrl,]
probe_groups$expected_channel <- "Green"
colnames(probe_groups) = c("address","control_probe","hypothetical_color","name","expected_channel")

# Merge with your intensity data
ctrlData <- ctrlData %>%
  left_join(probe_groups, by = "address")


# Everything is measured in the green channel and should be low
# As target removal should leave low residue signal

# Pivot to wide format (one row per address-sample pair)
ctrlDatal <- ctrlData %>%
  pivot_wider(
    names_from = channel,      # Column to split (Red/Green)
    values_from = value,       # Values to fill into new columns
    names_glue = "{tolower(channel)}"  # Force lowercase column names
  )

ctrlDatal <- data.frame(ctrlDatal)

sample_independent_controls_target_removal <- ggplot(ctrlDatal, aes(x = green, y = red, color = name)) +
  # Add diagonal x=y line first (so it's behind points)
  geom_abline(intercept = 0, slope = 1, color = "gray50", linetype = "dashed") +
  # Plot points
  geom_point(size = 1, alpha = 0.3) +
  # Make axes equal and set limits
  coord_equal(xlim = c(5,17), ylim = c(5, 17)) +
  labs(
    title = "Target removal controls",
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

sample_independent_controls_target_removal_by_array <- ggplot(ctrlData, aes(x = value, y = array, color = expected_channel)) + # Map color here for box outlines
  geom_boxplot(show.legend = FALSE) + # Use geom_boxplot and hide the legend for color
  scale_x_continuous(limits = c(5, 17)) +                 # Set the x-axis range from 5 to 17
  scale_color_manual(
    values = c("Green" = "darkgreen", "Red" = "darkred") # Assign specific colors for Green and Red channels
    # No 'name' argument needed here as the legend is hidden
  ) +
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


# Define bad samples by mean - 2 * sd
finite_values <- ctrlData$value[is.finite(ctrlData$value)]
mean_of_finite <- mean(finite_values, na.rm = TRUE)
sd_of_finite <- sd(finite_values, na.rm = TRUE)
bad_samples <- ctrlData[ctrlData$value > (mean_of_finite + 2 * sd_of_finite), ]$sample

# Write output
write.csv(
  na.omit(bad_samples),
  file = opt$flag_output,
  row.names = FALSE,
  col.names = FALSE,
)


pdf(opt$pdf_output, width = 8, height = 11)
grid.draw(ggplotGrob(sample_independent_controls_target_removal))
grid.newpage()
grid.draw(ggplotGrob(sample_independent_controls_target_removal_by_array))
grid.newpage()
grid.draw(ggplotGrob(sample_independent_controls_target_removal_by_slide))
dev.off()
