# Negative Control Probes: Design and Purpose
# ==========================================

#' Negative Control Probe Analysis
#' 
#' Evaluates system background noise using negative control probes in methylation arrays.

## Probe Design
# -------------
# - Sequences: Randomly permuted sequences designed NOT to hybridize to DNA template
# - Special importance in methylation studies due to reduced sequence complexity 
#   after bisulfite conversion
# - 600 negative control probes on array (Industry standard for Illumina platforms)

## Purpose and Function
# --------------------
# - Define system background signal (baseline noise floor)
# - GenomeStudio uses these to:
#   * Calculate average background signal 
#   * Determine standard deviation of background
#   * Establish detection limits for methylation probes
# - Must be evaluated in BOTH channels (Green/Red)

## Performance Metrics
# -------------------
# | Purpose Name       | Count | Green Channel | Red Channel | Expected Value |
# |--------------------|-------|---------------|-------------|-----------------|
# | Negative Average   | 600   |  (+)          |  (+)        | Background      |

## Usage
# ------
# Rscript 02_sample_dependent_controls_negative.R \
#   -r RGset.RData \
#   -p sample_dependent_controls_negative_overview.pdf

# Silent library loading and error handling
suppressPackageStartupMessages({
  tryCatch({
    library(minfi)
    library(tidyverse)
    library(reshape2)
    library(dplyr) 
    library(ggplot2)
    library(optparse)
  }, error = function(e) {
    stop("Package loading failed: ", e$message, call. = FALSE)
  })
})


option_list <- list(
  make_option(c("-r", "--rgset"), type="character", default=NULL,
              help="Path to RGset RData file"),
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


ctrl = "NEGATIVE"
red <- getRed(RGset)
green <- getGreen(RGset)
ctrlAddress <- getControlAddress(RGset, controlType = ctrl)
ctrlData <-  rbind(
            cbind(channel = "Red", melt(red[ctrlAddress, ], varnames = c("address", "sample"))),
            cbind(channel = "Green", melt(green[ctrlAddress, ], varnames = c("address", "sample"))))
ctrlData$value <- log2(ctrlData$value)
rm(red, green, ctrlAddress)

ctrlData <- ctrlData %>%
  group_by(sample, channel) %>%
  summarise(
    average = mean(value, na.rm = TRUE),
    std = sd(value, na.rm = TRUE),
    .groups = 'drop'
  )

ctrlData <- data.frame(ctrlData)

pdf(opt$pdf_output, width = 8, height = 11) 

ggplot(ctrlData, aes(x = channel, y = average, fill = channel)) +
  geom_boxplot(outlier.size = 0.5, alpha = 0.7) +
  scale_fill_manual(values = c("Red" = "red", "Green" = "green")) +
  labs(title = "Average intensity of 600 negative control probes per channel",
       x = "",
       y = "Intensity") +
  theme_minimal() +
  theme(legend.position = "none",
  plot.title = element_text(hjust = 0.5))

dev.off()
