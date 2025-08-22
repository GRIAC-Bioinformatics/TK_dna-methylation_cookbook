#!/usr/bin/env Rscript
# =============================================================================
# 02_meta_data_overview.R
# Author: Vartika Bisht
# Created: 13-Aug-2025
# Description:
# Generates comprehensive visual summaries of Illumina methylation array
# metadata including sample distributions across slides, arrays, plates,
# and wells. Provides quality control overviews for large-scale studies.
# Arguments:
#    --rgset           → Path to the RGChannelSet object file (.RData)
#    --pca             → Comma-separated list of variables for PCA (e.g., var1,var2,var3)
#    --stackbarplot    → Variables for stacked bar plots in the format grouping_variable:plotting_variable (e.g., var1:var2)
#    --platecol        → Column name for plate information (default: "Plate")
#    --slidecol        → Column name for slide information (default: "Slide")
#    --output          → Output PDF file path
# Usage:
#   Rscript 02_meta_data_overview.R --rgset <FILE.RData> \
#       --pca "Sample_ID,Plate" --stackbarplot "Slide:Plate" --output <output.pdf>
# Notes:
#   - This script requires a pre-processed RGChannelSet object as input.
#   - It helps identify potential batch effects or technical biases related to
#     sample processing and layout.
# =============================================================================

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
  make_option(c("-r", "--rgset"), type="character", default=NULL,
              help="Path to RGset.RData file", metavar="FILE"),
  make_option(c("-p", "--pca"), type="character", default=NULL,
              help="Variables for PCA (comma-separated, e.g., -p var1,var2,var3)", metavar="VAR1,VAR2,..."),
  make_option(c("-b", "--stackbarplot"), type="character", default=NULL,
              help="Variables for stacked bar plots comparisons in the format grouping_variable:plotting_variable (comma-separated, e.g., -b var1:var2,var2:var3,var3:var4)", metavar="VAR1:VAR2,VAR2:VAR3,..."),
  make_option(c("-l", "--platecol"), type="character", default="Plate",
              help="Column name for plate information", metavar="STRING"),
  make_option(c("-s", "--slidecol"), type="character", default="Slide",
              help="Column name for slide information", metavar="STRING"),
  make_option(c("-o", "--output"), type="character", default=NULL,
              help="Output PDF file path [default=%default]", metavar="FILE")
)

# Parse arguments
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

required_args <- c("rgset", "pca", "stackbarplot", "platecol", "slidecol", "output")

for (arg_name in required_args) {
  if (is.null(opt[[arg_name]])) {
    # If a required argument is NULL (not provided), print help and stop
    print_help(opt_parser)
    stop(paste("Error: Required argument --", arg_name, " is missing.", sep = ""))
  }
}

# Check if the rgset file exists
if (!file.exists(opt$rgset)) {
  stop(paste("Error: rgset file '", opt$rgset, "' not found.", sep = ""))
}

message("Loading RGset data from: ", opt$rgset)

# Load the RGset data
load(opt$rgset)

# Colnames
available_col_names <- colnames(RGset@colData)

# Validate --platecol
if (!(opt$platecol %in% available_col_names)) {
  stop(paste("Error: Specified plate column '", opt$platecol, "' not found in RGset colData. Available columns: ", 
             paste(available_col_names, collapse = ", "), sep = ""))
}

# Validate --slidecol
if (!(opt$slidecol %in% available_col_names)) {
  stop(paste("Error: Specified slide column '", opt$slidecol, "' not found in RGset colData. Available columns: ", 
             paste(available_col_names, collapse = ", "), sep = ""))
}

# Validate --pca format and variables
pca_vars <- unlist(strsplit(opt$pca, ","))
if (length(pca_vars) == 0) {
  stop("Error: --pca argument cannot be empty. Please provide comma-separated variables.")
}

invalid_pca_vars <- pca_vars[!(pca_vars %in% available_col_names)]
if (length(invalid_pca_vars) > 0) {
  stop(paste("Error: Invalid PCA variable(s) found: '", paste(invalid_pca_vars, collapse = ", "), 
             "'. These are not in RGset colData. Available columns: ", paste(available_col_names, collapse = ", "), sep = ""))
}

# Validate --stackbarplot format and variables
stack_barplot_pairs_str <- unlist(strsplit(opt$stackbarplot, ","))
if (length(stack_barplot_pairs_str) == 0) {
  stop("Error: --stackbarplot argument cannot be empty. Please provide comma-separated pairs (e.g., var1:var2).")
}

for (pair_str in stack_barplot_pairs_str) {
  parts <- unlist(strsplit(pair_str, ":"))
  
  if (length(parts) != 2) {
    stop(paste("Error: Invalid --stackbarplot format for '", pair_str, "'. Expected 'grouping_variable:plotting_variable'.", sep = ""))
  }
  
  grouping_var <- parts[1]
  plotting_var <- parts[2]
  
  if (!(grouping_var %in% available_col_names)) {
    stop(paste("Error: Grouping variable '", grouping_var, "' in --stackbarplot pair '", pair_str, 
               "' not found in RGset colData. Available columns: ", paste(available_col_names, collapse = ", "), sep = ""))
  }
  
  if (!(plotting_var %in% available_col_names)) {
    stop(paste("Error: Plotting variable '", plotting_var, "' in --stackbarplot pair '", pair_str, 
               "' not found in RGset colData. Available columns: ", paste(available_col_names, collapse = ", "), sep = ""))
  }
}

# --- All Validations Passed ---
message("All command-line arguments and data validations passed successfully. Proceeding with analysis.")

message("Creating output PDF: ", opt$output)

# Start PDF device
pdf(opt$output, 
    width = 12, 
    height = 8, 
    paper = "a4r")

## Enhanced function to create stacked bar plots
create_stacked_barplot <- function(data, title, x_label = "", y_label = "Proportion") {
  tryCatch({
    # Convert data to long format for ggplot
    long_data <- as.data.frame(data) %>%
      tibble::rownames_to_column(var = "Category") %>%
      tidyr::gather(key = "Group", value = "Count", -Category)
    
    long_data_proportions <- long_data %>%
      group_by(Category) %>% # Group by Category to calculate sum for each category
      mutate(Proportion = Count / sum(Count)) %>% # Calculate proportion
      ungroup()

    p <- ggplot(data.frame(long_data_proportions), aes(x = Category, y = Proportion, fill = Group)) +
          geom_bar(stat = "identity", position = "fill") + 
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

  # Calculate variance explained by each principal component
  pca_variance <- pca_result$sdev^2
  explained_variance_ratio <- pca_variance / sum(pca_variance) * 100

  # Create a data frame for the scree plot
  scree_data <- data.frame(
    PC = 1:10,
    Variance_Explained = explained_variance_ratio[1:10]
  )

  # Generate the elbow (scree) plot
  print(ggplot(scree_data, aes(x = PC, y = Variance_Explained)) +
          geom_bar(stat = "identity", fill = "steelblue", alpha = 0.7) + # Bar plot for individual variance
          geom_line(color = "red", size = 1) + # Line plot connecting points
          geom_point(color = "red", size = 3) + # Points on the line
          labs(
            title = "Scree Plot: Percentage of Variance Explained by PCs", # Main title
            x = "Principal Component",                                   # X-axis label
            y = "Percentage of Variance Explained (%)"                   # Y-axis label
          ) +
          theme_minimal() + # Use a minimal theme
          theme(
            plot.title = element_text(hjust = 0.5, face = "bold", size = 16), # Center and style title
            axis.title = element_text(size = 14),                             # Style axis titles
            axis.text = element_text(size = 12)                               # Style axis text
          ) +
          scale_x_continuous(breaks = scree_data$PC))

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
    # Split data by the grouping variable (e.g., Plate, Slide)
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
for (var in strsplit(opt$stackbarplot,",")[[1]]) {
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

