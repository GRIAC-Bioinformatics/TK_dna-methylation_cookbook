# Rscript: 02_sample_independent_controls_extension.R  -r RGset.RData -c controls_probes.csv -p sample_independent_controls_extension_overview.pdf -f Flagged_extension.csv
# Author: Vartika Bisht
# Date: 22 May 2025
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
# Extension controls test the extension efficiency of A, T, C and G from a hairpin probe, i.e.  probe contains its own complementary strand eliminating the need of DNA -- sample independent .
# A and T marked by DNP ( red channel ) and C and G marked by Biotin ( green channel )


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
              help="Output file for extension failures"),
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

rm(red, green, ctrlAddress)

probe_groups <- control_probes[control_probes$V2 == ctrl,]
probe_groups$V4 <- substr(probe_groups$V4,12,12)
probe_groups$V5 <- ifelse(probe_groups$V4 %in% c("A","T"),"Red","Green")
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
  )

summarized_df <- data.frame(ctrlData %>%
  group_by(slide, name, expected_channel) %>%
  summarise(value = median(value), .groups = 'drop')) %>%
  mutate(slide = as.factor(slide))
summarized_df <- summarized_df[is.finite(summarized_df$value),]

sample_independent_controls_extension_by_slide <- ggplot(summarized_df, aes(y = slide, x = value)) +
  geom_point(aes(color = expected_channel), size = 3) + # Map 'expected_channel' to color
  scale_x_continuous(limits = c(6, 17)) +                # Set x-axis limits from 6 to 15
  scale_color_manual(
    values = c("Green" = "darkgreen", "Red" = "darkred"), # Assign specific colors
    name = "Channel"                                     # Legend title for color
  ) +
  labs(
    title = "Extension controls by slide",
    x = "Median Intensity", # Update x-axis label
    y = "Slide ID"          # Update y-axis label
  ) +
  theme_minimal() + # A clean, minimal theme
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold") # Center and bold the title
  )

# Define bad samples by mean - 2 * sd
finite_values <- ctrlData$value[is.finite(ctrlData$value)]
mean_of_finite <- mean(finite_values, na.rm = TRUE)
sd_of_finite <- sd(finite_values, na.rm = TRUE)
bad_samples <- ctrlData[ctrlData$value < (mean_of_finite - 2 * sd_of_finite), ]$sample

# Write output
write.csv(
  na.omit(bad_samples),
  file = opt$flag_output,
  row.names = FALSE,
  col.names = FALSE,
)


pdf(opt$pdf_output, width = 8, height = 11)
grid.draw(ggplotGrob(sample_independent_controls_extension))
grid.newpage()
grid.draw(ggplotGrob(sample_independent_controls_extension_by_array))
grid.newpage()
grid.draw(ggplotGrob(sample_independent_controls_extension_by_slide))
dev.off()


