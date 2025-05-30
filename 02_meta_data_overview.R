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
# Rscript ${SCRIPTDIR}/02_meta_data_overview.R -i RGset.RData -o RGset_metadata_summary.pdf -p "Sample_Group,Sample_Plate,Sample_Well,Gender,Array,Slide" -l "Sample_Plate" -s "Slide" -b "Sample_Plate:Sample_Group,Sample_Plate:Gender,Slide:Sample_Group,Slide:Gender"
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
  make_option(c("-b", "--stackbarplot"), type="character", default=NULL,
              help="Variables for stacked bar plots comparisons in the format grouping_variable:plotting_variable (comma-separated, e.g., -b var1:var2,var2:var3,var3:var4)", metavar="VAR1:VAR2,VAR2:VAR3,..."),
  make_option(c("-l", "--platecol"), type="character", default="Sample_Plate",
              help="Column name for plate information", metavar="STRING"),
  make_option(c("-s", "--slidecol"), type="character", default="Slide",
              help="Column name for slide information", metavar="STRING"),
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

if (!is.null(opt$stackbarplot)) {
  stackbar_vars <- unlist(strsplit(opt$stackbarplot, ","))
} else {
  stackbar_vars <- NULL
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
            legend.position = "right") +
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
    paste("Object Size:", format(object.size(RGset), units = "GB")),
    paste("Unique Samples:", length(unique(RGset@colData$Sample_Name)))
  )

  # Add PCA variables if provided
  if (!is.null(pca_vars) && length(pca_vars) > 0) {
    summary_text <- c(summary_text, 
                      "PCA Variables:", 
                      paste(" -", pca_vars, collapse = "\n"))
    # Add sample plate and well information
    for (var in pca_vars) {
      summary_text <- c(summary_text,
                        paste("Unique", var, ":", 
                              length(unique(RGset@colData[[var]])) ) )
    }
  } else {
    summary_text <- c(summary_text, "No PCA variables provided.")
  }
  
  
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
        theme(legend.position = "right") +
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

plot_data_and_create_barplot <- function(rg_set, group_var, plot_var) {
  tryCatch({
    # Split data by the grouping variable (e.g., Sample_Plate, Slide)
    split_results <- lapply(split(rg_set@colData, rg_set@colData[[group_var]]),
                            function(x) table(x[[plot_var]]))
    
    # Convert list of tables to a list of data frames
    df_list <- lapply(names(split_results), function(group_val) {
      df <- as.data.frame(split_results[[group_val]])
      colnames(df) <- c(plot_var, group_val)
      df
    })
    
    # Merge all data frames into a single combined data frame
    combined_df <- Reduce(function(x, y) merge(x, y, by = plot_var, all = TRUE), df_list)
    rownames(combined_df) <- combined_df[[plot_var]]
    combined_df <- combined_df[,-1, drop = FALSE] # Remove the original variable column
    combined_df[is.na(combined_df)] <- 0 # Replace NA with 0
    
    # Transpose the data for plotting (rows become plates/slides, columns become categories)
    plot_data <- as.data.frame(t(combined_df))
    
    # Create the stacked bar plot
    create_stacked_barplot(plot_data, 
                          title = paste0(plot_var, " distribution by ", group_var),
                          x_label = group_var,
                          y_label = "")
  }, error = function(e) {
    message(paste0("Error creating ", plot_var, " by ", group_var, " plot: "), e$message)
  })
}

# Generate plots grouped by plate 
for (var in stackbar_vars) {
  # x axis variable is the grouping variable
  group_var <- strsplit(var, ":")[[1]][1]
  # This is the variable to be plotted on y axis
  plot_var <- strsplit(var, ":")[[1]][2]
  if (length(group_var) != 1 || length(plot_var) != 1) {
    stop("Invalid format for stackbarplot variable. Use 'var1:var2' format.")
  }
  message(paste("Creating stacked bar plot for", group_var, "by", plot_var))
  # Call the function to create the stacked bar plot
  plot_data_and_create_barplot(RGset, group_var , plot_var)
}

# Close the PDF device
dev.off()

message("Successfully created: ", opt$output)

