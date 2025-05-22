# Rscript: 02_sample_independent_controls_staining.R
# Author: Vartika Bisht
# Date: 22 May 2025

# --- Purpose ---
# This script processes sample-independent staining controls to assess
# the efficiency of staining in both the Red and Green channels.

# --- Staining Controls Probe Design ---
# These controls evaluate the effectiveness of the staining process.
# DNP (Dinitrophenol) immobilized beads are evaluated in the Red channel,
# while Biotin immobilized beads are evaluated in the Green channel.
# There are 4 probes in total: 2 for each color (Red/Green).

# | Staining Type  | Replicate | Evaluate Green (GRN) | Evaluate RED (RED) | Signal Level  |
# | :------------- | :-------- | :------------------- | :----------------- | :------------ |
# | DNP (High)     | 1         | -                    | +                  | High          |
# | DNP (Bgnd)     | 1         | -                    | +                  | Background    |
# | Biotin (Med)   | 1         | +                    | -                  | High          |
# | Biotin (Bgnd)  | 1         | +                    | -                  | Background    |

# --- Sample Independent Controls (STAINING) ---
# These are specific probes used for staining quality assessment:
#   - 21630339, STAINING, -99, DNP(20K), XXX
#   - 27630314, STAINING, Red, DNP (High),
#   - 43603326, STAINING, Purple, DNP (Bkg),
#   - 41666334, STAINING, Green, Biotin (High),
#   - 24669308, STAINING, -99, Biotin(5K), XXX
#   - 34648333, STAINING, Blue, Biotin (Bkg),

# --- Rationale for Analysis ---
# Staining is a critical step that occurs directly on the slide during hybridization.
# Therefore, it's crucial to assess staining efficiency for each individual **slide**
# and also for each **array** on the slide, as this reflects the hybridization quality.

# --- Input Parameters ---
# -r / --RGset_file: Path to the RData file containing the RGset object.
# -c / --controls_probes_file: Path to a CSV file detailing the control probe IDs

# --- Output Files ---
# -p / --output_pdf: Path and filename for the PDF report containing plots.
# -f / --flagged_csv: Path and filename for a CSV file listing any flagged
#                     samples.

# --- Usage ---
# Rscript 02_sample_independent_controls_staining.R \
#   -r RGset.RData \
#   -c controls_probes.csv \
#   -p sample_independent_controls_staining_overview.pdf \
#   -f Flagged_staining.csv

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
  make_option(c("-f", "--flag_output"), type="character", default=NULL,
              help="Output file for staining failures"),
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

probe_groups <- control_probes[control_probes$V2 == ctrl,]
probe_groups <- probe_groups[probe_groups$V3 != -99,]
probe_groups$expected_channel <- ifelse(grepl("DNP", probe_groups$V4), "Red", "Green")
colnames(probe_groups) = c("address","control_probe","hypothetical_color","name","expected_channel")

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
axis_max <- max(c(ctrlDatal$green, ctrlDatal$red), na.rm = TRUE)
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

summarized_df <- data.frame(ctrlData %>%
  group_by(slide, name) %>%
  summarise(value = median(value), .groups = 'drop'))
summarized_df <- summarized_df[is.finite(summarized_df$value),]
summarized_df <- summarized_df %>%
  mutate(
    type = ifelse(grepl("Biotin", name), "Biotin", "DNP"),
    level = ifelse(grepl("Bkg", name), "Bkg", "High")
  ) %>%
  mutate(slide = as.factor(slide))

sample_independent_controls_staining_by_slide <- ggplot(summarized_df, aes(y = slide, x = value, color = type, alpha = level)) + # Swap x and y aesthetics
  geom_point(size = 3) + 
  scale_color_manual(
    values = c("Biotin" = "darkgreen", "DNP" = "darkred"), # Define specific colors for Biotin and DNP
    name = "Probe Type"
  ) +
  scale_alpha_manual(
    values = c("High" = 1, "Bkg" = 0.6), # Define alpha (transparency) for High (darker) and Bkg (lighter but visible)
    name = "Level"
  ) +
  labs(
    title = "Staining controls by slide",
    x = "Median Intensity", # Update x-axis label
    y = "Slide ID"    # Update y-axis label
  ) +
  theme_minimal() + # A clean, minimal theme
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold")
  )


bad_samples = data.frame("V1"=union(ctrlData[ctrlData$name == "DNP (High)" & ctrlData$value < median(ctrlData[ctrlData$name == "DNP (Bkg)",]$value),]$sample,ctrlData[ctrlData$name == "Biotin (High)" & ctrlData$value < median(ctrlData[ctrlData$name == "Biotin (Bkg)",]$value),]$sample))


# Write output
write.csv(
  na.omit(bad_samples),
  file = opt$flag_output,
  row.names = FALSE,
  col.names = FALSE,
)


pdf(opt$pdf_output, width = 8, height = 11)
grid.draw(ggplotGrob(sample_independent_controls_staining))
grid.newpage()
grid.draw(ggplotGrob(sample_independent_controls_staining_by_array))
grid.newpage()
grid.draw(ggplotGrob(sample_independent_controls_staining_by_slide))
dev.off()


