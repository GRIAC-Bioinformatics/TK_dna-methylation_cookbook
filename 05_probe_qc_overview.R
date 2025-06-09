# Author: Vartika Bisht
# Date: 5 June 2025
#
# This R script is designed to perform quality control and summarization of DNA methylation probe data.
# It identifies and flags probes based on various criteria such as cross-reactivity,
# high detection p-values, location on sex chromosomes, high intensity, low bead count,
# and presence of SNPs. The script then generates an overview of these flagged probes,
# outputting a summary CSV file and a PDF plot with relevant visualizations.
#
# Rscript 05_probe_qc_overview.R --cross_reactive cross_reactive_flagged_probes.csv --detectionP Flagged_probe_detectionP.csv --high_intensity high_intensity_flagged_probes.csv --low_bead_count low_beadcount_flagged_probes.csv --snp_containing snp_containing_flagged_probes.csv --sex_chromosome 0 --output_csv_name probe_qc_overview.csv --output_plot_name probe_qc_overview.pdf
#
#
# Input:
#   --cross_reactive [FILE]: CSV file containing a list of cross-reactive probe IDs.
#   --detectionP [FILE]: CSV file containing a list of probe IDs with high detection p-values.
#   --high_intensity [FILE]: CSV file containing a list of probe IDs with high intensity.
#   --low_bead_count [FILE]: CSV file containing a list of probe IDs with low bead count.
#   --snp_containing [FILE]: CSV file containing a list of probe IDs with SNPs at SBE and CpG sites.
#   --sex_chromosome [NUM]: Numeric flag (0 or 1) to indicate whether to exclude (0) or include (1) sex chromosome probes.
#   --additional [FILE] (Optional): Path to an optional CSV file(s) for additional flagged probes.
#
# Output:
#   --output_csv_name [FILE]: Name for the output CSV file summarizing all flagged probes and their overlap.
#   --output_plot_name [FILE]: Name for the output PDF file containing an UpSet plot visualizing probe intersections.

suppressPackageStartupMessages({
  tryCatch({
library(optparse)
library(ggplot2)
library(minfi)
library(UpSetR)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
  }, error = function(e) {
    stop("Package loading failed: ", e$message, call. = FALSE)
  })
})

option_list <- list(
  make_option(c("-C", "--cross_reactive"), type = "character", default = NULL,
              help = "File containing list of cross reactive probes", metavar = "FILE"),
  make_option(c("-P", "--detectionP"), type = "character", default = NULL,
              help = "File containing list of probes with a high detectionP value", metavar = "FILE"),
  make_option(c("-S", "--sex_chromosome"), type = "numeric", default = 0,
              help = "1 or 0, whether to include(1) or exclude(0) sex chromosome probes", metavar = "NUM"),
  make_option(c("-H", "--high_intensity"), type = "character", default = NULL,
              help = "File containing list of probes with high intensity", metavar = "FILE"),
  make_option(c("-L", "--low_bead_count"), type = "character", default = NULL,
              help = "File containing list of probes with low bead count", metavar = "FILE"),
  make_option(c("-N", "--snp_containing"), type = "character", default = NULL,
              help = "File containing list of probes containing SNPs at SBE and CpG sites", metavar = "FILE"),
  make_option(c("-a", "--additional"), type = "character", default = NULL,
              help = "Path(s) to an optional CSV file(s) for Additional probes flagged."),
  make_option(c("-f", "--min_flag_overlap"), type = "character", default = NULL,
              help = "Minimum number of overlapping flags required to filter a sample (required)."),
  make_option(c("-o", "--output_csv_name"), type = "character", default = NULL,
              help = "Name for the output CSV file (required)."),
  make_option(c("-p", "--output_plot_name"), type = "character", default = NULL,
              help = "Name for the output PDF plot file containing all plots (required).")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Validate input arguments
if (is.null(opt$cross_reactive) || !file.exists(opt$cross_reactive)) {
  stop("Please provide the path to the cross reactive probes file using -C or --cross_reactive option.")
}
if (is.null(opt$detectionP) || !file.exists(opt$detectionP)) {
  stop("Please provide the path to the detectionP probes file using -P or --detectionP option.")
}
if (is.null(opt$high_intensity) || !file.exists(opt$high_intensity)) {
  stop("Please provide the path to the high intensity probes file using -H or --high_intensity option.")
}
if (is.null(opt$low_bead_count) || !file.exists(opt$low_bead_count)) {
  stop("Please provide the path to the low bead count probes file using -L or --low_bead_count option.")
}
if (is.null(opt$snp_containing) || !file.exists(opt$snp_containing)) {
  stop("Please provide the path to the SNP containing probes file using -N or --snp_containing option.")
}
if (is.null(opt$output_csv_name)) {
  stop("Please provide the path to save the output CSV file using -o or --output_csv_name option.")
}
if (is.null(opt$output_plot_name)) {
  stop("Please provide the path to save the PDF file using -p or --output_plot_name option.")
}

# Load the flagged probes data
cross_reactive_probes <- as.character(read.csv(opt$cross_reactive, header = TRUE, stringsAsFactors = FALSE)[,1])
detectionP_probes <- as.character(read.csv(opt$detectionP, header = TRUE, stringsAsFactors = FALSE)[,1])
high_intensity <- as.character(read.csv(opt$high_intensity, header = TRUE, stringsAsFactors = FALSE)[,1])
low_bead_count <- as.character(read.csv(opt$low_bead_count, header = TRUE, stringsAsFactors = FALSE)[,1])
snp_containing <- as.character(read.csv(opt$snp_containing, header = TRUE, stringsAsFactors = FALSE)[,1])

if(is.null(opt$additional)){
    additional <- c()
} else {
    additional <- as.character(read.csv(opt$additional, header = TRUE, stringsAsFactors = FALSE)[,1])
}

sex_chromosome <- as.numeric(opt$sex_chromosome)
if(sex_chromosome!=1){
    message("Sex chromosome to be excluded")
    anno <- minfi::getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
    sex_chromosome <- anno$Name[anno$chr %in% c("chrX", "chrY")]
} else {
    message("Sex chromosome to be included")
    sex_chromosome <- c()
}

message("Read all flagged probes..")

# --- Combine all your lists into a single named list ---
all_probes_raw_input <- list(
  "Cross_Reactive" = cross_reactive_probes,
  "Detection_P" = detectionP_probes,
  "High_Intensity" = high_intensity,
  "Low_Bead_Count" = low_bead_count, # This is our intended empty list for testing
  "SNP_Containing" = snp_containing,
  "Sex_Chromosome" = sex_chromosome,
  "Additional" = additional # Your new list
)

filtered_probes_list <- list()

# Iterate through the raw input list
for (list_name in names(all_probes_raw_input)) {
  current_list <- all_probes_raw_input[[list_name]]
  if (length(current_list) == 0) {
    message(paste0("Warning: The list '", list_name, "' is empty and will be excluded from the UpSet plot."))
  } else {
    filtered_probes_list[[list_name]] <- current_list
  }
}

# Check if there are any lists left to plot
if (length(filtered_probes_list) == 0) {
  stop("Error: All lists are empty. Cannot generate an UpSet plot.")
}

message("Generating UpSet plot..")

# --- Generate the UpSet plot ---
# The nsets parameter should now be dynamic based on the filtered list
num_sets_to_plot <- length(filtered_probes_list)

pdf(opt$output_plot_name, height = 8 , width = 11)
print(upset(
  fromList(filtered_probes_list),
  nsets = num_sets_to_plot,
  nintersects = NA, # Plot all possible intersections
  order.by = "freq", # Order intersections by frequency
  decreasing = TRUE, # Show highest frequency first
  sets.bar.color = "#56B4E9", # Customize set size bar color
  mainbar.y.label = "Intersection Size",
  sets.x.label = "Set Size",
  text.scale = c(1.3, 1.3, 1, 1, 2, 0.75), # Adjust text sizes for better readability
  point.size = 3, # Size of dots in the matrix
  line.size = 1 # Thickness of lines connecting dots
))
dev.off()

# --- Create and save the binary matrix as a CSV ---
# Get all unique probe names across all filtered lists
all_unique_probes <- unique(unlist(filtered_probes_list))

# Initialize an empty data frame with probes as row names and lists as column names
binary_matrix_df <- data.frame(
  row.names = all_unique_probes,
  matrix(0,
         nrow = length(all_unique_probes),
         ncol = length(filtered_probes_list),
         dimnames = list(all_unique_probes, names(filtered_probes_list)))
)

# Fill the data frame with 1s where a probe is present in a list
for (list_name in names(filtered_probes_list)) {
  probes_in_list <- filtered_probes_list[[list_name]]
  binary_matrix_df[probes_in_list, list_name] <- 1
}

# This sums the '1's across each row for all set columns
binary_matrix_df$Total_Occurrence <- rowSums(binary_matrix_df)

# Save the binary matrix to a CSV file
write.csv(binary_matrix_df, opt$output_csv_name, row.names = TRUE)
message("Flagged probes from filters saved to : ", opt$output_csv_name)

