# Rscript: 02_sample_independent_controls_hybridization.R
# Author: Vartika Bisht
# Date: 23 May 2025

# --- Purpose ---
# This script is designed to assess the overall performance of the entire assay
# by analyzing hybridization controls, which utilize synthetic targets instead of amplified DNA.

# --- Hybridization Controls Details ---
# The hybridization controls utilize synthetic targets that perfectly complement the sequences on the array,
# allowing the probe to extend on these synthetic targets as a template.
# These synthetic targets are present in the Hybridization Buffer (RA1) at three distinct concentration levels:
# high-concentration (5 pM), medium-concentration (1 pM), and low-concentration (0.2 pM).
# All bead type IDs associated with these controls are expected to result in signals with varying intensities,
# directly corresponding to the concentrations of the initial synthetic targets.
# It is important to note that the performance of hybridization controls should *only* be monitored in the green channel.

# --- Hybridization Controls Design ---
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

# --- Input Parameters ---
# -r / --RGset_file: Path to the RData file containing the hybridizationRGset objecthybridization.
# -c / --controls_probes_file: Path to a CSV file that details the hybridizationcontrol probe IDshybridization

# --- Output Files ---
# -p / --output_pdf: Path and filename for the generated hybridizationPDF reporthybridization. 
#                    (e.g., sample_independent_controls_hybridization_overview.pdf)
# -f / --flagged_csv: Path and filename for a hybridizationCSV filehybridization that will list any flagged
#                     samples or arrays. 

# --- Usage ---
# To execute this script from your R environment or terminal, use the following command:
# Rscript 02_sample_independent_controls_hybridization.R \
#   -r RGset.RData \
#   -c controls_probes.csv \
#   -p sample_independent_controls_hybridization_overview.pdf \
#   -f Flagged_hybridization.csv


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
              help="Output file for hybridization failures"),
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

ctrl = "HYBRIDIZATION"
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
probe_groups$V4 <- substr(probe_groups$V4,6,nchar(probe_groups$V4)-1)
probe_groups$expected_channel <- "Green"
colnames(probe_groups) = c("address","control_probe","hypothetical_color","name","expected_channel")

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

sample_independent_controls_hybridization <- ggplot(ctrlDatal, aes(x = green, y = red, color = name)) +
  # Add diagonal x=y line first (so it's behind points)
  geom_abline(intercept = 0, slope = 1, color = "gray50", linetype = "dashed") +
  # Plot points
  geom_point(size = 1, alpha = 0.3) +
  # Make axes equal and set limits
  coord_equal(xlim = c(5,17), ylim = c(5, 17)) +
  scale_color_manual(
    values = c(
      "High" = "#006400",    # DarkGreen
      "Medium" = "#32CD32",  # LimeGreen
      "Low" = "#90EE90"      # LightGreen
    ),
    name = "Probe Type" # Keep the legend title
  ) +
  labs(
    title = "Hybridization controls",
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

summarized_df <- data.frame(ctrlData %>%
  group_by(slide, name) %>%
  summarise(value = median(value), .groups = 'drop'))  %>%
  mutate(slide = as.factor(slide))
summarized_df <- summarized_df[is.finite(summarized_df$value),]

sample_independent_controls_hybridization_by_slide <- ggplot(summarized_df, aes(y = slide, x = value, color = name)) + # Swap x and y aesthetics
  geom_point(size = 3) + 
   scale_color_manual(
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
    title = "Target removal controls by slide",
    x = "Median Intensity", # Update x-axis label
    y = "Slide ID"    # Update y-axis label
  ) +
  theme_minimal() + # A clean, minimal theme
  scale_x_continuous(limits = c(5, 17)) +  
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold")
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
  file = opt$flag_output,
  row.names = FALSE,
  col.names = FALSE,
)


pdf(opt$pdf_output, width = 8, height = 11)
grid.draw(ggplotGrob(sample_independent_controls_hybridization))
grid.newpage()
grid.draw(ggplotGrob(sample_independent_controls_hybridization_by_array))
grid.newpage()
grid.draw(ggplotGrob(sample_independent_controls_hybridization_by_slide))
dev.off()
