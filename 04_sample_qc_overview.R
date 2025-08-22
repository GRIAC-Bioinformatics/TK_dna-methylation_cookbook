#!/usr/bin/env Rscript
# =============================================================================
# 04_sample_qc_overview.R
# Author: Vartika Bisht
# Created: 13-Aug-2025
# Description:
# This script generates a comprehensive overview of all samples flagged during
# various upstream quality control (QC) steps. It combines flagging information
# from multiple QC scripts into a single summary report, providing plots and a
# combined data file to identify samples with significant quality issues.
# Arguments:
#    --rgset                 → Path to the RGChannelSet object file (.RData)
#    --platform              → DNA methylation platform (e.g., EPIC, 450K)
#    --manifestkey           → Path to the manifest key file (default: data/manifest.annotation.key.csv)
#    --assembly              → Genome assembly version
#    --pca_vars              → Comma-separated list of metadata variables for PCA and pie charts
#    --pca_rdata             → Output RData file for PCA results
#    --flaggedcombined       → Output CSV file with all combined flag information
#    --flaggedlist           → Input CSV file listing the QC flags and paths to their respective flagged sample files
#    --pdf_output            → Output PDF file containing all QC overview plots
#    --min_flag_overlap      → Minimum number of flags a sample must have to be included in PCA (default: 2)
# Usage:
#   Rscript 04_sample_qc_overview.R --rgset <FILE.RData> --pca_vars <VAR1,VAR2> \
#     --flaggedlist <flag_paths.csv> --flaggedcombined <combined_flags.csv> \
#     --pca_rdata <pca_results.RData> --pdf_output <qc_overview.pdf>
# Notes:
#   - This script integrates flags from various quality control scripts, including:
#     - Sample-Dependent Controls: Bisulfite Conversion, Specificity, Nonpolymorphic, Negative
#     - Sample-Independent Controls: Staining, Extension, Target Removal, Hybridization
#     - Additional Controls: Sex Mismatch, Detection P-value, Low Intensity
#   - The combined flag information is used to generate visualizations, such as PCA, to help
#     identify samples with systematic technical issues or batch effects.
#   - The min_flag_overlap parameter allows for filtering out samples with only minor
#     QC issues, focusing the analysis on samples with more severe quality problems.
# =============================================================================


suppressPackageStartupMessages({
  tryCatch({
library(optparse)
library(ggplot2)
library(minfi)
library(UpSetR)
library(grid) 
library(methods) 
library(reshape2) 
library(dplyr)
library(ggrepel)
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
  make_option(c("-v", "--pca_vars"), type = "character", default = NULL,
              help = "Comma-separated list of variables from RGset@colData to use for PCA coloring and pie charts (required)."),
  make_option(c("-d","--pca_rdata"), type = "character", default = NULL,
              help = "Name for the output RData file containing PCA results (required)."),
  make_option(c("-c", "--flaggedcombined"), type = "character", default = NULL,
              help = "Name for the output CSV file with all flag information combined (required)."),
  make_option(c("-f", "--flaggedlist"), type = "character", default = NULL,
              help = "CSV file with list of QC flag and the path to the CSV file with the list of samples that failed that QC step (required)."),
  make_option(c("-o", "--pdf_output"), type = "character", default = NULL,
              help = "Name for the output PDF plot file containing all plots (required)."),
  make_option(c("-m", "--min_flag_overlap"), type = "numeric", default = 2,
              help = "Minimum number of flags a sample must have to be included in PCA and pie charts (default: 2).")
)

# Parse command line options
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)


# VALIDATE REQUIRED ARGUMENTS
required_args <- c("rgset", "platform", "manifestkey", "assembly" , "pca_vars" , "flaggedcombined" , "flaggedlist" , "pdf_output", "min_flag_overlap","pca_rdata")

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


# Load the RGset data
message("Loading RGset data from: ", opt$rgset)
load(opt$rgset)


message("Load CSV file with the list of QC flags and the samples that failed each QC step from ", opt$flaggedlist, " ...")
flaggedlist <- tryCatch({
  read.csv(opt$flaggedlist, header = FALSE)
}, error = function(e) {
  stop("Failed to read flagged list file: ", e$message)
})

# The first column of flagged.list is the QC flag name, and the second column is the path to the CSV file with the list of samples that failed that QC step
# Make a list of list with the list name as the QC flag name and the list of samples that failed that QC step
flagged <- list()
for (i in 1:nrow(flaggedlist)) {
  flag <- flaggedlist$V1[i]
  path <- flaggedlist$V2[i]
  message("Reading flagged samples for flag: ", flag)
  flag_samples <- tryCatch({
    read.csv(path, header = TRUE)[[1]]
  }, error = function(e) {
    stop("Failed to read flagged samples file for flag '", flag, "': ", e$message)
  })

  if(length(flag_samples) == 0) {
    message("No samples found for flag: ", flag)
  } else {
    message("Found ", length(flag_samples), " samples for flag: ", flag)
    flagged[[flag]] <- flag_samples
  }
}

message("Prosessing flags : ", paste(names(flagged), collapse = " , "))

# Create a binary data.frame with samples as rows and flags as columns
# Each cell will be TRUE if the sample belongs to that flag, FALSE otherwise
all_flagged_samples <- unique(unlist(flagged))
all_flags <- as.character(names(flagged))

# Initialize empty data.frame of FALSE
flagged_df <- as.data.frame(
  matrix(0, nrow = length(all_flagged_samples), ncol = length(all_flags)),
  row.names = all_flagged_samples
)
colnames(flagged_df) <- all_flags

# Fill in TRUE where sample belongs to a flag
for (flag in all_flags) {
  message("Reading flagged samples for flag: ", flag)
  samplelist <- as.character(flagged[[flag]])
  flagged_df[samplelist, flag] <- 1
}

# Add 'Total Occurrences' column (based on combined categories)
flagged_df$Total_Occurrences <- rowSums(flagged_df)
flagged_df$Sample <- rownames(flagged_df)

# Sort the table by 'Total_Occurrences' in descending order
flagged_df <- flagged_df[order(flagged_df$Total_Occurrences, decreasing = TRUE),c("Sample", "Total_Occurrences", all_flags)]

write.csv(flagged_df, file = opt$flaggedcombined, row.names = FALSE) 
message(paste("Sample QC overview saved to:", opt$flaggedcombined))


# Plotting Section 
pdf(opt$pdf_output, width = 12, height = 8) # Open PDF device for all plots

# Apply `as.numeric` to convert these columns from character/logical to numeric (0 or 1)
flagged_df[all_flags] <- lapply(flagged_df[all_flags], as.numeric)

message("Generate UpSet Plot...")

# Generate the UpSet Plot
# The UpSet plot is designed to show the intersections of multiple sets (your QC flags).
# - The top bars represent the size of each *intersection* (i.e., how many samples fall into that specific combination of flags).
# - The dots and lines matrix below indicates which flags are part of each intersection.
# - The bars on the left show the total number of samples affected by *each individual flag* (set size).
tryCatch({
  if (nrow(flagged_df) > 0 && length(all_flags) > 0) {
    # Specify which columns in your data frame are the 'sets' for the UpSet plot
    sets_to_plot <- all_flags

    # Create the UpSet plot
    print(upset(
      flagged_df,               # Your data frame containing 0s and 1s for flags
      sets = sets_to_plot,        # The names of the columns to treat as sets
      main.bar.color = "steelblue", # Color for the main intersection bars (top)
      sets.bar.color = "darkgreen", # Color for the individual set size bars (left)
      matrix.color = "darkblue",    # Color for the dots in the intersection matrix
      point.size = 2.5,             # Size of the dots in the matrix
      line.size = 1,                # Thickness of the lines connecting dots in the matrix
      text.scale = c(1.3, 1.3, 1, 1, 1.5, 1.2), # Adjust font sizes for different plot elements
      order.by = "freq",            # Order intersections by their frequency (most common first)
      empty.intersections = "on",   # Include combinations that have zero samples
      mainbar.y.label = "Intersection Size (Number of Samples)", # Y-axis label for top bars
      sets.x.label = "Set Size (Total Samples per Flag)" # X-axis label for left bars
    ))

    message("UpSet Plot generated successfully.")

  } else {
    warning("No samples or flag columns found for UpSet plot.")
    message("Please ensure 'flagged_df' has rows and relevant flag columns (e.g., 'Flagged_...csv').")
  }
}, error = function(e) {
  message("Error generating UpSet Plot: ", e$message)
  message("Please check your data and ensure 'UpSetR' is installed and loaded correctly.")
})


# 2. & 3. Load RGset and Filter Samples
message("Filtering samples for PCA and Pie Charts...")

tryCatch({
  # Identify samples flagged at least twice based on flagged_df
  # These are the samples to EXCLUDE from PCA
  min_flag_overlap <- as.numeric(opt$min_flag_overlap)
  flagged_samples <- flagged_df$Sample[flagged_df$Total_Occurrences >= min_flag_overlap]
  occurrences_count <- flagged_df$Total_Occurrences[flagged_df$Total_Occurrences >= min_flag_overlap]
  names(occurrences_count) <- flagged_samples

  message(paste("Found", length(flagged_samples), "samples flagged at least ",min_flag_overlap," times (will be highlighted in the PCA)."))

  pca_vars <- strsplit(opt$pca_vars,",")[[1]]
  if (length(pca_vars) == 0) {
    warning("No valid PCA variables remaining after checking RGset@colData for samples not flagged >= ",min_flag_overlap," times. Skipping PCA.")
    grid.newpage()
    grid.text("No valid PCA variables for plotting samples not flagged >= ",min_flag_overlap," times.", gp = gpar(col = "red"))
  } else {
    # Principal Components Analysis (PCA) plot (color flagged samples)
    message("Generating PCA Plots for samples flagged at least ", min_flag_overlap, " times...")
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
      
      save(pca_result, file = opt$pca_rdata)
      message("PCA completed successfully. PCA RData file saved at: ", opt$pca_rdata)

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
      message("Error creating PCA plot for samples not flagged >= ",min_flag_overlap," times: ", e$message)
      grid.newpage()
      grid.text(paste("Error displaying PCA plot for samples not flagged >= ",min_flag_overlap," times:", e$message),
                gp = gpar(col = "red"))
    })
  }}) 


  tryCatch({
    # Pie Charts for flagged samples (>= min_flag_overlap times)
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
message(paste("All plots saved to:", opt$pdf_output))