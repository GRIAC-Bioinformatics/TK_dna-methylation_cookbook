#!/usr/bin/env Rscript
# =============================================================================
# 03_sample_sex_mismatch.R
# Author: Vartika Bisht
# Created: 13-Aug-2025
# Description:
# This script performs a sex mismatch quality control check. It predicts the
# sex of each sample based on copy number information from the X and Y
# chromosomes and compares this prediction to the sex provided in the metadata.
# This serves as a key metric for detecting potential sample mix-ups.
# Arguments:
#    --grset                → Path to the GenomicRatioSet object file (.RData)
#    --flagged              → Output CSV file path for flagged samples with a sex mismatch
#    --pdf_output           → Output PDF file path for QC plots
# Usage:
#   Rscript 03_sample_sex_mismatch.R --grset <FILE.RData> \
#     --flagged <sex_mismatch_failures.csv> --pdf_output <output.pdf>
# Notes:
#   - The getSex function from minfi is used to predict sex based on the
#     median total intensity (copy number) of probes on the X and Y chromosomes.
#   - Predicted Sex: The difference between median Y and X chromosome intensities
#     (yMed - xMed) is a bimodal distribution.
#       - Females (XX): Have a higher median X and very low median Y intensity,
#         resulting in a negative yMed - xMed value.
#       - Males (XY): Have a lower median X and a detectable median Y intensity,
#         resulting in a positive yMed - xMed value.
#   - A sex mismatch occurs when the predicted sex does not align with the sex
#     specified in the metadata. This is a strong indicator of a sample swap.
#   - The script visualizes this by plotting a histogram of the yMed - xMed
#     values, with bars colored to show whether a sample's predicted sex matches
#     its metadata. A red bar indicates a mismatch.
# Aryee MJ, Jaffe AE, Corrada-Bravo H, Ladd-Acosta C, Feinberg AP, Hansen KD, Irizarry RA. Minfi: a flexible and 
# comprehensive Bioconductor package for the analysis of Infinium DNA methylation microarrays. Bioinformatics. 2014 
# May 15;30(10):1363-9. doi: 10.1093/bioinformatics/btu049. Epub 2014 Jan 28. PMID: 24478339; PMCID: PMC4016708.
# =============================================================================

# Silent library loading and error handling
suppressPackageStartupMessages({
  tryCatch({
  library(minfi)
  library(ggplot2)
  library(dplyr)
  library(optparse)
  }, error = function(e) {
    stop("Package loading failed: ", e$message, call. = FALSE)
  })
})

# Set up command line options
option_list <- list(
  make_option(c("-g", "--grset"), type="character", default=NULL,
              help="Path to GRset file", metavar="FILE"),
  make_option(c("-f", "--flagged"), type="character", default=getwd(),
              help="CSV file with flagged samples where sex mismatch", metavar="DIR"),
  make_option(c("-p", "--pdf_output"), type="character", default=getwd(),
              help="PDF file with all plots", metavar="DIR")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

# VALIDATE REQUIRED ARGUMENTS
required_args <- c("grset", "flagged" , "pdf_output")

for (arg in required_args) {
  if (is.null(opt[[arg]])) {
    stop(paste("Error: Required argument --", arg, " is missing. Use --help for more information.", sep = ""))
  }
}


message("Checking sex mismatch...")

load(opt$grset)
# Observed sex is taken from the column in the colData of the GRset object
message("Recording the sex as defined in the metadata sheet.")
observedSex <- grSet@colData[["Gender"]]
names(observedSex) <- rownames(grSet@colData)

message("Pridicting Sex ..")
message("Predicted Sex is determined using the copy number values of X and Y chromosomes")
message("Get copy number values which are defined as the sum of the methylation and unmethylation channel.")
CN <- getCN(grSet)

# Calculate median values for X and Y chromosomes
xIndex <- which(seqnames(grSet) == "chrX")
yIndex <- which(seqnames(grSet) == "chrY")
xMed <- colMedians(CN, rows = xIndex, na.rm = TRUE, useNames = TRUE)
yMed <- colMedians(CN, rows = yIndex, na.rm = TRUE, useNames = TRUE)

# Calculate the difference between Y and X medians
message("Calculate the difference in the copy number variation between X and Y chromosomes.")
message("The difference should seperate into 2 clusters. One which is closer to 0, other closer to -4.")
message("Highly negative associates to 2 X chromosome , i.e. female and close to 0 is for males.")
dd <- yMed - xMed

# Remove samples with NaNs which can occur if there are no probes on the X or Y chromosome
# This can happen if the sample are controls.
message("Removing following samples as they have NaN values from the difference...")
message(length(dd) , " rows before removing NA..")
dd <- dd[!is.na(dd)]
dd <- data.frame((dd))
colnames(dd) <- "CopyNumberDifference"
dd <- cbind(dd, observedSex = observedSex[rownames(dd)])
message(nrow(dd) , " rows after removing NA..")

# remove samples with no sex information
dd <- dd[dd$observedSex != "-", ]

# Find local maxima to determine centers for clustering
message("Finding local maxima to determine centers for clustering...")

local_max_indices <- c()  

for (i in seq_along(density(dd$CopyNumberDifference)$y)) {
  if (i == 1 || i == length(density(dd$CopyNumberDifference)$y)) {
    next
  }
  if (density(dd$CopyNumberDifference)$y[i] > density(dd$CopyNumberDifference)$y[i - 1] && 
      density(dd$CopyNumberDifference)$y[i] > density(dd$CopyNumberDifference)$y[i + 1]) {
    local_max_indices <- c(local_max_indices, i)
  }
}

female_center <- min(density(dd$CopyNumberDifference)$x[local_max_indices])
male_center <- max(density(dd$CopyNumberDifference)$x[local_max_indices])

message("Local maxima indices found: ", paste(density(dd$CopyNumberDifference)$x[local_max_indices], collapse = ", "))
message("Female center: ", female_center)
message("Male center: ", male_center)

# Check variation and run clustering
if (length(unique(dd$CopyNumberDifference)) > 2) {
  message("Clustering samples based on copy number values...")
  message("Number of samples with valid copy number values: ", length(dd$CopyNumberDifference))
  message("Centers for clustering will be set to ", female_center, " for females and ", male_center, " for males.")
  # Clusetring using k-means
  k <- kmeans(dd$CopyNumberDifference, c(female_center,male_center), iter.max = 10, nstart = 1, algorithm = "Hartigan-Wong")
  # Determine which cluster - min for females (1) - max for males (2)
  predictedSex <- k$cluster
  predictedSex <- ifelse(k$cluster == 1 , "F","M")

    
  dd <- cbind(dd, predictedSex)
  dd$MatchStatus <- ifelse(dd$observedSex == dd$predictedSex, "Match", "Mismatch")

  message("Generating plots and saving to PDF file: ", opt$pdf_output)
  pdf(opt$pdf_output, width = 6, height = 6)
  message("Plotting density and histogram of copy number differences...")
  print(ggplot(dd, aes(x = CopyNumberDifference)) +
    geom_density(fill = "#1E90FF", color = "#333333", alpha = 0.8, linewidth = 0.8) + 
    labs(title = "Distribution of Sex Chromosome\nCopy Number Differences", 
        x = expression(bold("Copy Number Difference (Y - X Chromosome)")),
        y = expression(bold("Density"))) + 
    theme_minimal() +
    theme(
      panel.grid.major = element_line(color = "#E0E0E0", linewidth = 0.5), 
      panel.grid.minor = element_blank(), 
      panel.border = element_blank(), 
      legend.position = "bottom"
    ) +
    geom_vline(xintercept=female_center, linetype="dashed", color="black", linewidth=1) +
    geom_vline(xintercept=male_center, linetype="dashed", color="black", linewidth=1))

  message("Plotting histogram of copy number differences with match status...")
  print(ggplot(dd, aes(x = CopyNumberDifference, fill = MatchStatus)) +
    geom_histogram( alpha = 0.7, binwidth=0.05) + 
    labs(title = "Distribution of Sex Chromosome\nCopy Number Differences", 
        x = expression(bold("Copy Number Difference (Y - X Chromosome)")),
        y = expression(bold("Frequency"))) + 
    theme_minimal() +
    theme(
      panel.grid.major = element_line(color = "#E0E0E0", linewidth = 0.5), 
      panel.grid.minor = element_blank(), 
      panel.border = element_blank(), 
      legend.position = "bottom"
    ) +
    scale_fill_manual(values = c("Match" = "#90EE90", "Mismatch" = "#FF0000")) + 
    geom_vline(xintercept=female_center, linetype="dashed", color="black", linewidth=1) +
    geom_vline(xintercept=male_center, linetype="dashed", color="black", linewidth=1))


  beta <- getBeta(grSet)
    
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
  
  message("Plotting PCA with sex mismatch samples highlighted")
  plot_data_pca$MatchStatus <- "Match"
  plot_data_pca[rownames(dd[dd$MatchStatus == "Mismatch",]), "MatchStatus"] <- "Mismatch"
  print(ggplot(plot_data_pca, aes(x = PC1, y = PC2), color = "gray50") +
    geom_point(alpha = 0.7) +
    geom_point(data=plot_data_pca[plot_data_pca$MatchStatus == "Mismatch",], aes(x = PC1, y = PC2), color = "red") +
    theme(legend.position = "right") +
    labs( title = "PCA with all samples",
          subtitle = "Samples with sex mismatch are highlighted",
          x = paste0("PC1 (", pca_var_per[1], "%)"),
          y = paste0("PC2 (", pca_var_per[2], "%)")) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold")))

  dev.off()
  message("PDF saved to ", opt$pdf_output)
  write.csv(data.frame("sexmismatch" = rownames(dd[dd$MatchStatus == "Mismatch",])), file = opt$flagged, row.names = FALSE)
  message("Samples with sample sex mismatch are saved to ", opt$flagged)
} else {
  warning("Not enough variation to estimate sex..")
  print(dd)
  write.csv(data.frame(Sample = character()), file = opt$flagged, row.names = FALSE)
}
