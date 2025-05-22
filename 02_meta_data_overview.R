#!/usr/bin/env Rscript
# meta.data.overview.R - Methylation Array Metadata Visualization Script
#
# Author: Vartika Bisht
# Date:   20 May 2025
#
# Purpose: Generate comprehensive visual summaries of Illumina methylation array
#          metadata including sample distributions across slides, arrays, plates,
#          and wells. Provides quality control overviews for large-scale studies.
#
# Input:   RGset.RData file containing RGChannelSet object
# Output:  PDF report with visualizations of sample metadata distributions
#
# Parameters:
#   -i, --input    Path to input RGset.RData file (required)
#   -p, --pca      Variables for PCA (space-separated, e.g., -p var1 var2 var3)
#   -o, --output   Output PDF path [default: ./RGset_metadata_summary.pdf]
# Usage: Rscript meta.data.overview.R -i input.RData -o report.pdf -p Sample_Group Sample_Plate Sample_Well Gender Array Slide
#
# Hierarchy Documentation:
#   Slide:       Physical chip identifier (e.g., 6264488065)
#   Array:       Position on slide (R01C01 = Row 01 Column 01)
#   Sample_Plate: Processing plate identifier (e.g., MeDALL1)
#   Sample_Well: Well position on plate (e.g., A01, B01)
#
# Note: Each sample is uniquely identified by the combination of:
#       Slide + Array + Plate + Well identifiers
#       (e.g., 6264488065 + R01C01 + MeDALL1 + A01)

# Silent library loading and error handling
suppressPackageStartupMessages({
  tryCatch({
    library(minfi)
    library(ggplot2)
    library(tidyr)
    library(dplyr)
    library(gridExtra)
    library(grid)
    library(optparse)
  }, error = function(e) {
    stop("Package loading failed: ", e$message, call. = FALSE)
  })
})


# Set up command line options
option_list <- list(
  make_option(c("-i", "--input"), type="character", default=NULL,
              help="Path to RGset.RData file", metavar="FILE"),
  make_option(c("-p", "--pca"), type="character", default=NULL,
              help="Variables for PCA (comma-separated, e.g., -p var1,var2,var3)", metavar="VAR1,VAR2,..."),
  make_option(c("-o", "--output"), type="character", default=NULL,
              help="Output PDF file path [default=%default]", metavar="FILE")
)

# Parse arguments
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

# Check for required input
if (is.null(opt$input)) {
  print_help(opt_parser)
  stop("Input file must be specified (-i or --input)", call.=FALSE)
}

# Verify input file exists
if (!file.exists(opt$input)) {
  stop(paste("Input file not found:", opt$input), call.=FALSE)
}

message("Loading RGset data from: ", opt$input)

# Load the RGset data
load(opt$input)

# Verify the loaded object
if (!exists("RGset") || !inherits(RGset, "RGChannelSet")) {
  stop("The loaded object is not an RGChannelSet or is not named 'RGset'", call.=FALSE)
}

if (!is.null(opt$pca)) {
  pca_vars <- unlist(strsplit(opt$pca, ","))
} else {
  stop("No PCA variables provided. Use -p var1,var2 ...")
}

message("Creating output PDF: ", opt$output)

# Start PDF device
pdf(opt$output, 
    width = 12, 
    height = 8, 
    paper = "a4r")

## Enhanced function to create stacked bar plots
create_stacked_barplot <- function(data, title, x_label = "", y_label = "Count") {
  tryCatch({
    # Convert data to long format for ggplot
    long_data <- as.data.frame(data) %>%
      tibble::rownames_to_column(var = "Category") %>%
      tidyr::gather(key = "Group", value = "Count", -Category)
    
    # Create the plot
    p <- ggplot(long_data, aes(x = Category, y = Count, fill = Group)) +
      geom_bar(stat = "identity", position = "stack") +
      labs(title = title, x = x_label, y = y_label) +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            plot.title = element_text(hjust = 0.5, face = "bold"),
            legend.position = "bottom") +
      scale_fill_brewer(palette = "Dark2")
    
    # Print the plot
    grid.newpage()
    grid.draw(ggplotGrob(p))
    
  }, error = function(e) {
    message("Error creating bar plot: ", e$message)
    grid.newpage()
    grid.text(paste("Error displaying plot:", title), 
              gp = gpar(col = "red"))
  })
}

## Title Page (unchanged)
tryCatch({
  grid.newpage()
  grid.text("RGset Metadata Summary", 
            gp = gpar(fontsize = 24, fontface = "bold"), 
            x = 0.5, y = 0.7)
  grid.text(paste("Generated on", Sys.Date()), 
            gp = gpar(fontsize = 18), 
            x = 0.5, y = 0.6)
  grid.text("Comprehensive Methylation Array Analysis", 
            gp = gpar(fontsize = 16), 
            x = 0.5, y = 0.5)
}, error = function(e) {
  message("Error creating title page: ", e$message)
})

## Text Summary Page (unchanged)
tryCatch({
  grid.newpage()
  summary_text <- c(
    paste("Dataset Dimensions:", paste(dim(RGset), collapse = " x ")),
    paste("Object Size:", format(object.size(RGset), units = "GB")),
    paste("Unique Samples:", length(unique(RGset@colData$Sample_Name))),
    paste("Unique Slides:", length(unique(RGset@colData$Slide))),
    paste("Unique Arrays:", length(unique(RGset@colData$Array))),
    paste("Unique Sample Plates:", length(unique(RGset@colData$Sample_Plate))),
    paste("Unique Sample Groups:", length(unique(RGset@colData$Sample_Group))),
    paste("Unique Genders:", length(unique(RGset@colData$Gender))),
    "Sample Group Distribution:",
    capture.output(print(table(RGset@colData$Sample_Group))),
    "",
    "Gender Distribution:",
    capture.output(print(table(RGset@colData$Gender))),
    ""
  )
  
  # Split long text into multiple pages if needed
  text_lines <- unlist(strsplit(paste(summary_text, collapse = "\n"), "\n"))
  lines_per_page <- 50
  text_pages <- split(text_lines, ceiling(seq_along(text_lines)/lines_per_page))
  
  for (i in seq_along(text_pages)) {
    if (i > 1) grid.newpage()
    page_title <- if (length(text_pages) > 1) {
      paste("Data Summary (Page", i, "of", length(text_pages), ")")
    } else {
      "Data Summary"
    }
    grid.text(page_title, x = 0.05, y = 0.97, just = c("left", "top"),
              gp = gpar(fontsize = 12, fontface = "bold"))
    grid.text(paste(text_pages[[i]], collapse = "\n"), 
              x = 0.05, y = 0.93, 
              just = c("left", "top"),
              gp = gpar(fontsize = 10, fontfamily = "mono"))
  }
}, error = function(e) {
  message("Error creating summary page: ", e$message)
})


## Principal Components Analysis (PCA) plot
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
  plot_data <- data.frame(
    PC1 = pca_result$x[, 1],
    PC2 = pca_result$x[, 2],
    as.data.frame(RGset@colData)[, pca_vars]
  )
  
  # Modified plotting function
  plot_pca_by_variable <- function(data, variable) {
    n_unique <- length(unique(data[[variable]]))
    
    if (n_unique > 8) {
      # Plot text labels instead of points for high-cardinality variables
      ggplot(data, aes(x = PC1, y = PC2, label = .data[[variable]])) +
        geom_text(size = 3, alpha = 0.7, check_overlap = TRUE) +
        labs(title = paste("PCA -", variable),
             x = paste0("PC1 (", pca_var_per[1], "%)"),
             y = paste0("PC2 (", pca_var_per[2], "%)")) +
        theme_minimal() +
        theme(legend.position = "none")
    } else {
      # Regular plot with points and legend for low-cardinality variables
      ggplot(data, aes(x = PC1, y = PC2, color = .data[[variable]])) +
        geom_point(alpha = 0.7) +
        labs(title = paste("PCA colored by", variable),
             x = paste0("PC1 (", pca_var_per[1], "%)"),
             y = paste0("PC2 (", pca_var_per[2], "%)")) +
        theme_minimal() +
        theme(legend.position = "bottom") +
        guides(color = guide_legend(nrow = 2, byrow = TRUE))
    }
  }
  
  # Create and display one plot per page
  for (var in pca_vars) {
    print(plot_pca_by_variable(plot_data, var))
  }
  
}, error = function(e) {
  message("Error creating PCA plot: ", e$message)
  grid.newpage()
  grid.text(paste("Error displaying PCA plot:", e$message), 
            gp = gpar(col = "red"))
})

## NEW: Create Female vs Male per Plate plot
tryCatch({
  # Create count of females and males per plate
  gender_by_plate <- as.data.frame.matrix(
    table(RGset@colData$Sample_Plate, RGset@colData$Gender)
  )
  
  create_stacked_barplot(gender_by_plate, 
                        "Gender Distribution by Plate", 
                        x_label = "Sample Plate")
}, error = function(e) {
  message("Error creating Gender by Plate plot: ", e$message)
})

## Sample Group by Plate plot
tryCatch({
  split_results <- lapply(split(RGset@colData, RGset@colData$Sample_Plate), 
                         function(x) table(x$Sample_Group))
  df_list <- lapply(names(split_results), function(plate) {
    df <- as.data.frame(split_results[[plate]])
    colnames(df) <- c("Sample_Group", plate)
    df
  })
  combined_df <- Reduce(function(x, y) merge(x, y, by = "Sample_Group", all = TRUE), df_list)
  rownames(combined_df) <- combined_df$Sample_Group
  combined_df <- combined_df[,-1, drop = FALSE]
  combined_df[is.na(combined_df)] <- 0
  
  # Transpose for plotting
  plot_data <- as.data.frame(t(combined_df))
  
  create_stacked_barplot(plot_data, 
                        "Sample Group Distribution by Plate", 
                        x_label = "Sample Plate")
}, error = function(e) {
  message("Error creating Sample Group by Plate plot: ", e$message)
})

## Sample Group by Slide plot
tryCatch({
  split_results <- lapply(split(RGset@colData, RGset@colData$Slide), 
                         function(x) table(x$Sample_Group))
  df_list <- lapply(names(split_results), function(slide) {
    df <- as.data.frame(split_results[[slide]])
    colnames(df) <- c("Sample_Group", slide)
    df
  })
  combined_df <- Reduce(function(x, y) merge(x, y, by = "Sample_Group", all = TRUE), df_list)
  rownames(combined_df) <- combined_df$Sample_Group
  combined_df <- combined_df[,-1, drop = FALSE]
  combined_df[is.na(combined_df)] <- 0
  
  # Transpose for plotting
  plot_data <- as.data.frame(t(combined_df))
  
  create_stacked_barplot(plot_data, 
                        "Sample Group Distribution by Slide", 
                        x_label = "Slide")
}, error = function(e) {
  message("Error creating Sample Group by Slide plot: ", e$message)
})

## Gender by Slide plot
tryCatch({
  split_results <- lapply(split(RGset@colData, RGset@colData$Slide), 
                         function(x) table(x$Gender))
  df_list <- lapply(names(split_results), function(slide) {
    df <- as.data.frame(split_results[[slide]])
    colnames(df) <- c("Gender", slide)
    df
  })
  combined_df <- Reduce(function(x, y) merge(x, y, by = "Gender", all = TRUE), df_list)
  rownames(combined_df) <- combined_df$Gender
  combined_df <- combined_df[,-1, drop = FALSE]
  combined_df[is.na(combined_df)] <- 0
  
  # Transpose for plotting
  plot_data <- as.data.frame(t(combined_df))
  
  create_stacked_barplot(plot_data, 
                        "Gender Distribution by Slide", 
                        x_label = "Slide")
}, error = function(e) {
  message("Error creating Gender by Slide plot: ", e$message)
})

# Close the PDF device
dev.off()

message("Successfully created: ", opt$output)

