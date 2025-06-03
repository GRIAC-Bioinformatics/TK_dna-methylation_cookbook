# R Script for Sample QC and Upset Plot
# Author: Vartika Bisht
# Date: 24 May 2025

# Introduction:
# This script processes multiple CSV files, each containing a single row of sample names.
# It consolidates these names into a single data frame, indicating the presence of each
# sample in different experimental categories (flags). A 'Total Occurrences' column
# summarizes how many times each sample appeared across all categories.
# The script then generates an 'UpSet' plot to visualize the intersections of samples
# across these categories, providing a comprehensive overview of sample overlap.
# The processed data is saved as 'sample_QC_overview.csv'.

# Usage:
# To run this script from the command line, use Rscript followed by the script name
# and the required (and optional) arguments.
# Note: For flags that can take multiple CSV files, provide them as a space-separated
# string enclosed in double quotes.
# Example:
# Rscript 03_sample_qc_overview.R \
#     -b "Flagged_bisulfite_conversion_I.csv Flagged_bisulfite_conversion_II.csv" \
#     -n "Flagged_nonpolymeric.csv" \
#     -c "Flagged_specificity_I.csv Flagged_specificity_II.csv" \
#     -e "Flagged_extension.csv" \
#     -h "Flagged_hybridization.csv" \
#     -s "Flagged_staining.csv" \
#     -t "Flagged_target_removal.csv" \
#     -G "RGset.Rdata" \
#     -P "Sample_Group Sample_Plate Sample_Well Gender Array Slide" \
#     -o "sample_qc.csv" \
#     -p "sample_qc_analysis.pdf" \
#     -a "optional_additional_file.csv optional_another_file.csv"

# Input:
# -b, --bisulfite_flag: (Required) Path(s) to the CSV file(s) for Bisulfite Flag samples.
# -n, --nonpolymeric_flag: (Required) Path(s) to the CSV file(s) for Nonpolymeric Flag samples.
# -c, --specificity_flag: (Required) Path(s) to the CSV file(s) for Specificity Flag samples.
# -e, --extension_flag: (Required) Path(s) to the CSV file(s) for Extension Flag samples.
# -y, --hybridisation_flag: (Required) Path(s) to the CSV file(s) for Hybridisation Flag samples.
# -s, --staining_flag: (Required) Path(s) to the CSV file(s) for Staining Flag samples.
# -t, --target_removal_flag: (Required) Path(s) to the CSV file(s) for Target Removal Flag samples.
# -G, --rgset: (Required) Path to the RGChannelSet .Rdata file (e.g., 'RGset.Rdata').
# -P, --pca_vars: (Required) Space-separated list of variables from RGset@colData to use for PCA coloring and pie charts.
# -a, --additional: (Optional) Path(s) to additional CSV file(s) for samples.

# Output:
# 1. sample_QC_overview.csv: A CSV file containing a summary table.
#    - Rows: Unique sample names found across all input files.
#    - Columns: Binary flags (0 or 1) indicating presence in each category,
#               and a 'Total_Occurrences' column showing the sum of flags for each sample.
#    - Sorted: By 'Total_Occurrences' in descending order.
# 2. sample_QC_upset_plot.png: A PNG image file of the UpSet plot.
#    - Visualizes the intersections of samples across different categories.
#    - X-axis: Represents the combinations of categories (intersections).
#    - Y-axis: Shows the frequency (number of samples) in each intersection.

suppressPackageStartupMessages({
  tryCatch({
library(optparse)
library(ggplot2)
library(minfi)
library(grid) 
library(methods) 
library(reshape2) 
library(dplyr)
library(ggrepel)
  }, error = function(e) {
    stop("Package loading failed: ", e$message, call. = FALSE)
  })
})


# 1. Define command-line options using optparse
option_list <- list(
  make_option(c("-b", "--bisulfite_flag"), type = "character", default = NULL,
              help = "Path(s) to the CSV file(s) for Bisulfite Flag samples (required)."),
  make_option(c("-n", "--nonpolymeric_flag"), type = "character", default = NULL,
              help = "Path(s) to the CSV file(s) for Nonpolymeric Flag samples (required)."),
  make_option(c("-c", "--specificity_flag"), type = "character", default = NULL,
              help = "Path(s) to the CSV file(s) for Specificity Flag samples (required)."),
  make_option(c("-e", "--extension_flag"), type = "character", default = NULL,
              help = "Path(s) to the CSV file(s) for Extension Flag samples (required)."),
  make_option(c("-y", "--hybridisation_flag"), type = "character", default = NULL,
              help = "Path(s) to the CSV file(s) for Hybridisation Flag samples (required)."),
  make_option(c("-s", "--staining_flag"), type = "character", default = NULL,
              help = "Path(s) to the CSV file(s) for Staining Flag samples (required)."),
  make_option(c("-t", "--target_removal_flag"), type = "character", default = NULL,
              help = "Path(s) to the CSV file(s) for Target Removal Flag samples (required)."),
  make_option(c("-G", "--rgset"), type = "character", default = NULL,
              help = "Path to the RGChannelSet .Rdata file (e.g., 'RGset.Rdata') (required)."),
  make_option(c("-P", "--pca_vars"), type = "character", default = NULL,
              help = "Space-separated list of variables from RGset@colData to use for PCA coloring and pie charts (required)."),
  make_option(c("-a", "--additional"), type = "character", default = NULL,
              help = "Path(s) to an optional CSV file(s) for Additional samples."),
  make_option(c("-o", "--output_csv_name"), type = "character", default = NULL,
              help = "Name for the output CSV file (required)."),
  make_option(c("-p", "--output_plot_name"), type = "character", default = NULL,
              help = "Name for the output PDF plot file containing all plots (required).")
)


# Parse command-line arguments
# The crucial fix: add_help_option=FALSE to prevent default -h flag creation
opt_parser <- OptionParser(option_list = option_list, add_help_option = FALSE)
opt <- parse_args(opt_parser)

# Check if all required arguments are provided
required_flags <- c("bisulfite_flag", "nonpolymeric_flag", "specificity_flag",
                    "extension_flag", "hybridisation_flag", "staining_flag",
                    "target_removal_flag", "rgset", "pca_vars",
                    "output_csv_name", "output_plot_name")

for (flag in required_flags) {
  if (is.null(opt[[flag]])) {
    stop(paste("Error: The --", flag, " flag is required. Please provide a path or name.", sep = ""))
  }
}

# Parse PCA variables
pca_vars <- unlist(strsplit(opt$pca_vars, " "))
if (length(pca_vars) == 0) {
  stop("Error: No PCA variables provided with -P/--pca_vars.")
}


# Function to read samples from a CSV file and return a data frame with source file
read_samples_with_source <- function(file_path) {
  if (!file.exists(file_path)) {
    warning(paste("Warning: File not found:", file_path, ". Skipping this file."))
    return(data.frame(sample_name = character(0), exact_filename = character(0), stringsAsFactors = FALSE))
  }
  tryCatch({
    samples_df_raw <- read.csv(file_path, header = TRUE, stringsAsFactors = FALSE)
    samples <- as.character(samples_df_raw[,1])
    samples <- samples[!is.na(samples) & samples != ""]
    
    # Return a data frame with sample names and the exact source filename
    return(data.frame(sample_name = samples, exact_filename = basename(file_path), stringsAsFactors = FALSE))
  }, error = function(e) {
    warning(paste("Error reading file:", file_path, ". Error:", e$message, ". Skipping this file."))
    return(data.frame(sample_name = character(0), exact_filename = character(0), stringsAsFactors = FALSE))
  })
}

# Initialize a data frame to store all raw sample data
all_samples_df_raw <- data.frame(
  sample_name = character(0),
  category_name = character(0), # This will be the combined category (e.g., "Bisulfite")
  exact_filename = character(0), # This will be the specific file name (e.g., "Flagged_bisulfite_conversion_I.csv")
  stringsAsFactors = FALSE
)

# Process each flag's input files
process_flag_files <- function(flag_value, category) {
  if (!is.null(flag_value)) {
    # Split the input string by spaces to get individual file paths
    file_paths <- unlist(strsplit(flag_value, " "))
    for (f_path in file_paths) {
      current_samples <- read_samples_with_source(f_path)
      if (nrow(current_samples) > 0) {
        current_samples$category_name <- category # Assign the combined category name
        all_samples_df_raw <<- rbind(all_samples_df_raw, current_samples)
      }
    }
  }
}

# Define the mapping from flag to combined category name
flag_category_map <- list(
  bisulfite_flag = "Bisulfite",
  nonpolymeric_flag = "Nonpolymeric",
  specificity_flag = "Specificity",
  extension_flag = "Extension",
  hybridisation_flag = "Hybridisation",
  staining_flag = "Staining",
  target_removal_flag = "Target_Removal",
  additional = "Additional"
)

# Call process_flag_files for each category based on the map
for (flag_name in names(flag_category_map)) {
  process_flag_files(opt[[flag_name]], flag_category_map[[flag_name]])
}

# Get all unique sample names
all_unique_samples <- unique(all_samples_df_raw$sample_name)

# If no samples are found, exit gracefully
if (length(all_unique_samples) == 0) {
  stop("No samples found in any of the provided CSV files. Please check your input files.")
}

# Get all unique exact filenames
all_unique_filenames <- unique(all_samples_df_raw$exact_filename)

# Create the main sample QC data frame
sample_qc_df <- data.frame(Sample = all_unique_samples, stringsAsFactors = FALSE)

# Add binary columns for each individual file
for (file_name in all_unique_filenames) {
  samples_in_file <- as.integer(sample_qc_df$Sample %in% all_samples_df_raw[all_samples_df_raw$exact_filename == file_name,]$sample)
  sample_qc_df[[file_name]] <- samples_in_file
}

# Add 'Total Occurrences' column (based on combined categories)
sample_qc_df$Total_Occurrences <- rowSums(sample_qc_df[, all_unique_filenames, drop = FALSE])

# Sort the table by 'Total_Occurrences' in descending order
sample_qc_df <- sample_qc_df[order(sample_qc_df$Total_Occurrences, decreasing = TRUE), ]

# Save the sample QC overview table
output_csv_filename <- opt$output_csv_name
write.csv(sample_qc_df, file = output_csv_filename, row.names = FALSE) 
message(paste("Sample QC overview saved to:", output_csv_filename))

# 3.  Plotting Section 
output_plot_filename <- opt$output_plot_name
pdf(output_plot_filename, width = 12, height = 8) # Open PDF device for all plots



# 1. Bar Plot: Samples vs. Total Occurrences
message("Generating Bar Plot...")
tryCatch({
  if (nrow(sample_qc_df) > 0) {
    p_bar <- ggplot(sample_qc_df, aes(x = reorder(Sample, -Total_Occurrences), y = Total_Occurrences)) +
      geom_bar(stat = "identity", fill = "steelblue") +
      theme_minimal() +
      labs(title = "Total QC Flags per Sample",
           x = "Sample",
           y = "Number of Flags") +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 8),
            plot.title = element_text(hjust = 0.5, face = "bold"))
    print(p_bar)
  } else {
    warning("No samples found for bar plot.")
    grid.newpage()
    grid.text("No samples available to generate bar plot.", gp = gpar(col = "red"))
  }
}, error = function(e) {
  message("Error generating Bar Plot: ", e$message)
  grid.newpage()
  grid.text(paste("Error generating Bar Plot:", e$message), gp = gpar(col = "red"))
})


# 2. & 3. Load RGset and Filter Samples
message("Loading RGset and filtering samples for PCA and Pie Charts...")
RGset_path <- opt$rgset
RGset <- NULL

if (!file.exists(RGset_path)) {
  stop(paste("Error: RGChannelSet Rdata file not found at:", RGset_path))
}

tryCatch({
  load(RGset_path) # This should load an object named 'RGset' into the environment
  if (!exists("RGset") || !is(RGset, "RGChannelSet")) {
    stop("The loaded Rdata file does not contain an object named 'RGset' of class 'RGChannelSet'.")
  }
  message("RGChannelSet loaded successfully.")

  # Identify samples flagged at least twice based on sample_qc_df
  # These are the samples to EXCLUDE from PCA
  flagged_samples <- sample_qc_df$Sample[sample_qc_df$Total_Occurrences >= 2]
  occurrences_count <- sample_qc_df$Total_Occurrences[sample_qc_df$Total_Occurrences >= 2]
  names(occurrences_count) <- flagged_samples

  message(paste("Found", length(flagged_samples), "samples flagged at least twice (will be highlighted in the PCA)."))

  if (length(pca_vars) == 0) {
    warning("No valid PCA variables remaining after checking RGset@colData for samples not flagged >= 2 times. Skipping PCA.")
    grid.newpage()
    grid.text("No valid PCA variables for plotting samples not flagged >= 2 times.", gp = gpar(col = "red"))
  } else {
    # Principal Components Analysis (PCA) plot (color flagged samples)
    message("Generating PCA Plots for samples not flagged at least twice...")
    tryCatch({
      # Extract methylation values
      beta <- getBeta(RGset)
        
      # Handle missing/infinite values
      if (any(is.na(beta))) {
        message("NA values found - imputing with row medians")
        beta_imputed <- apply(beta, 1, function(x) {
          x[is.na(x)] <- median(x, na.rm = TRUE)
          return(x)
        })
        beta <- t(beta_imputed)
      }
      
      if (any(!is.finite(beta))) {
        message("Infinite values found - replacing with row max/min")
        beta_fixed <- apply(beta, 1, function(x) {
          x[x == Inf] <- max(x[is.finite(x)], na.rm = TRUE)
          x[x == -Inf] <- min(x[is.finite(x)], na.rm = TRUE)
          return(x)
        })
        beta <- t(beta_fixed)
      }
      
      # Transpose and perform PCA
      beta_t <- t(beta)
      pca_result <- prcomp(beta_t, scale. = TRUE, center = TRUE)

      # Calculate variance explained
      pca_var <- pca_result$sdev^2
      pca_var_per <- round(pca_var/sum(pca_var)*100, 1)

      # Prepare plot data
      plot_data_pca <- data.frame(
        PC1 = pca_result$x[, 1],
        PC2 = pca_result$x[, 2]
      )
      rownames(plot_data_pca) <- rownames(pca_result$x)

      plot_data_pca$Total_Occurrences <- 0
      plot_data_pca[flagged_samples, "Total_Occurrences"] <- occurrences_count[flagged_samples]
      plot_data_pca$Total_Occurrences <- as.factor(plot_data_pca$Total_Occurrences)

      pcaplot <- ggplot(plot_data_pca, aes(x = PC1, y = PC2), color = "gray50") +
            geom_point(alpha = 0.7) +
            geom_point(data=plot_data_pca[plot_data_pca$Total_Occurrences != 0,], aes(x = PC1, y = PC2, color = Total_Occurrences)) +
            theme(legend.position = "right") +
            labs( title = "PCA with all samples",
                  subtitle = "Flagged samples are highlighted",
                  x = paste0("PC1 (", pca_var_per[1], "%)"),
                  y = paste0("PC2 (", pca_var_per[2], "%)"),
                  color = "Total Occurrences") +
          theme_minimal() +
          theme(plot.title = element_text(hjust = 0.5, face = "bold"))

      print(pcaplot)

    }, error = function(e) {
      message("Error creating PCA plot for samples not flagged >= 2 times: ", e$message)
      grid.newpage()
      grid.text(paste("Error displaying PCA plot for samples not flagged >= 2 times:", e$message),
                gp = gpar(col = "red"))
    })
  }}) 


  tryCatch({
    # Pie Charts for flagged samples (>= 2 times)
    message("Generating Pie Charts for samples flagged at least twice...")
    # Check if RGset@colData contains the required PCA variables
    pca_vars <- pca_vars[pca_vars %in% colnames(RGset@colData)]
    if (length(pca_vars) == 0) {
      warning("No valid PCA variables found in RGset@colData for pie charts. Skipping Pie Charts.")
      grid.newpage()
      grid.text("No valid PCA variables for plotting pie charts.", gp = gpar(col = "red"))
      return()
    }
    # Check if RGset@colData contains the required PCA variables
    if(length(pca_vars[!pca_vars %in% colnames(RGset@colData)])> 0) {
      warning(paste("The following PCA variables are not found in RGset@colData for pie charts:",
                    paste(pca_vars[!pca_vars %in% colnames(RGset@colData)], collapse = ", "), ". Skipping them for Pie Charts."))
      grid.newpage()
      grid.text(paste("The following PCA variables are not found in RGset@colData for pie charts:",
                    paste(pca_vars[!pca_vars %in% colnames(RGset@colData)], collapse = ", "), ". Skipping them for Pie Charts."), gp = gpar(col = "red"))
      return()
    }
    
    # Filter RGset@colData for samples flagged at least twice
    rgset_pie_subset <- RGset@colData[flagged_samples,pca_vars]
    if (dim(rgset_pie_subset)[1] == 0 || dim(rgset_pie_subset)[2] == 0) {
        warning("Either no samples flagged at least twice found in RGset@colData for pie charts or no pca variables.")
        grid.newpage()
        grid.text("Either no samples flagged at least twice found in RGset@colData for pie charts or no pca variables.", gp = gpar(col = "red"))
      } else {
        for (var in pca_vars) {
          # Count occurrences of each category
          counts <- table(rgset_pie_subset[[var]])
          df_pie <- as.data.frame(counts)
          colnames(df_pie) <- c("Category", "Count")
          df_pie$Percentage <- round(df_pie$Count / sum(df_pie$Count) * 100, 1)

          # Create a custom label for the legend that includes Category and Percentage
          df_pie$LegendLabel <- paste0(df_pie$Category, " (", df_pie$Percentage, "%)")

          p_pie <- ggplot(df_pie, aes(x = "", y = Count, fill = LegendLabel)) + # Map LegendLabel to fill
            geom_bar(width = 1, stat = "identity") +
            coord_polar("y", start = 0) +
            labs(
              title = paste("Distribution of", var, "\nin samples flagged at least once"),
              fill = "Category (Percentage)" # Title for the legend
            ) +
            theme_void() + # Minimal theme for pie chart
            theme(
              plot.title = element_text(hjust = 0.5, face = "bold"),
              legend.position = "right", # Place the legend on the right
              legend.text = element_text(size = 10), # Adjust legend text size
              legend.title = element_text(face = "bold") # Make legend title bold
            )

          print(p_pie)
        }
      }

  }, error = function(e) {
    message("Error processing RGset: ", e$message)
    grid.newpage()
    grid.text(paste("Error loading or processing RGset:", e$message), gp = gpar(col = "red"))
  })

dev.off() # Close the PDF device
message(paste("All plots saved to:", output_plot_filename))